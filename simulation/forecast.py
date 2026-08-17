# simulation/forecast.py

"""
Prospective forecast: predictions issued before the outcome is known.

WHY THIS IS DIFFERENT FROM THE HINDCAST
---------------------------------------
Every number this project has produced so far was measured against data that
already existed. Even the held-out split (`validation.py`) chose where to cut
*after* seeing the whole record — an honest test, but a retrospective one.

A forecast committed to version control before the outcome resolves cannot be
tuned, cannot leak, and cannot have its split chosen to flatter it. The commit
is the timestamp. That makes this the only test here that is impossible to
game, and the only one that can genuinely falsify the model.

WHAT IS BEING PREDICTED
-----------------------
The model's last observation is 2024, so everything from 2025 onward is out of
sample. That splits the forecast into three bands with different evidential
value:

    2025        already elapsed, resolvable immediately by anyone with the
                data. The model has never seen it.
    2026        partly elapsed at issue time.
    2027-2030   genuinely prospective.

PREDICTING OUR OWN FAILURES
---------------------------
Each prediction carries the subsystem's measured out-of-sample skill, so the
forecast states in advance which of its own numbers it expects to be wrong.
Solar tested at R2 = +0.555 and should verify; oceanic tested at -2.705 and
should not. If that pattern holds, the validation framework is itself
validated. If oceanic verifies well, the held-out result was misleading and
`validation.py` needs re-examining.

Intervals come from the HELD-OUT test RMSE, not the in-sample fit. Using
in-sample error would understate them by roughly a factor of two on oceanic.

WHAT THIS MODEL CANNOT DO
-------------------------
CEED has no ENSO mode. A very strong El Nino was developing at issue time, and
none of its effect is in these numbers. The atmospheric prediction is
therefore expected to run COLD for 2026-2027 — a directional, falsifiable
statement recorded here rather than explained afterwards. See E9 in
legacy/README.md for why the oscillator that would have supplied it was
rejected as the wrong shape.
"""

import argparse
import json
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional

import numpy as np

from simulation.convergence_model import (
    SYSTEMS,
    AnthropogenicForcing,
    SystemParameters,
)
from simulation.hindcast import integrate_from
from simulation.unit_bridge import ALL_OBSERVED, SCALES
from simulation.validation import DEFAULT_TRAIN_END, evaluate_split

LAST_OBSERVED_YEAR = 2024
DEFAULT_HORIZON = 6

# Interval width in units of the held-out test RMSE.
INTERVAL_SIGMA = 2.0


@dataclass
class SkillPrior:
    """Measured out-of-sample skill for one subsystem, from validation.py."""
    system: str
    test_r2: float
    test_rmse: float
    units: str

    @property
    def expectation(self) -> str:
        """What we expect of this subsystem's forecast, stated in advance."""
        if self.test_r2 > 0.30:
            return "expect skill"
        if self.test_r2 > 0.0:
            return "marginal"
        return "expect NO skill (flat line beats it)"


@dataclass
class Prediction:
    """One pre-registered number."""
    system: str
    year: int
    value: float
    lower: float
    upper: float
    units: str
    test_r2: float
    expectation: str
    band: str          # 'resolvable now' | 'resolving' | 'prospective'


def skill_priors(train_end: int = DEFAULT_TRAIN_END) -> Dict[str, SkillPrior]:
    """Measure out-of-sample skill per subsystem, live.

    Computed rather than hardcoded so the forecast can never quote a stale
    skill figure for a model that has since changed.
    """
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    split = evaluate_split(params, "skill priors", train_end=train_end)
    return {
        name: SkillPrior(system=name,
                         test_r2=score.r_squared,
                         test_rmse=score.rmse,
                         units=score.units)
        for name, score in split.test.items()
    }


def _band(year: int, issue_year: int) -> str:
    if year < issue_year:
        return "resolvable now"
    if year == issue_year:
        return "resolving"
    return "prospective"


def run_forecast(start_year: int = LAST_OBSERVED_YEAR,
                 horizon: int = DEFAULT_HORIZON,
                 issue_year: int = 2026,
                 params: Optional[SystemParameters] = None,
                 priors: Optional[Dict[str, SkillPrior]] = None
                 ) -> List[Prediction]:
    """Integrate forward from the last observation and emit predictions.

    Args:
        start_year: Year to initialise from. Must have observations.
        horizon: Years to forecast beyond start_year.
        issue_year: Year the forecast is issued, for banding.
        params: Model parameters. Defaults to the shipped anthropogenic set.
        priors: Skill priors. Computed if omitted.

    Returns:
        One Prediction per subsystem per forecast year.

    Raises:
        ValueError: if start_year has no observations to initialise from.
    """
    if start_year not in ALL_OBSERVED['solar']:
        raise ValueError(
            f"no observations at {start_year}; cannot initialise a forecast")

    if params is None:
        params = SystemParameters(anthropogenic=AnthropogenicForcing())
    if priors is None:
        priors = skill_priors()

    t, solution = integrate_from(start_year, horizon, params)

    out: List[Prediction] = []
    for year in range(start_year + 1, start_year + horizon + 1):
        t_target = float(year - start_year)
        for i, name in enumerate(SYSTEMS):
            E = float(np.interp(t_target, t, solution[:, i]))
            value = float(SCALES[name].energy_to_obs(E))
            prior = priors[name]
            half = INTERVAL_SIGMA * prior.test_rmse
            out.append(Prediction(
                system=name, year=year, value=value,
                lower=value - half, upper=value + half,
                units=prior.units,
                test_r2=prior.test_r2,
                expectation=prior.expectation,
                band=_band(year, issue_year),
            ))
    return out


def print_forecast(predictions: List[Prediction], issue_year: int = 2026):
    """Print the forecast grouped by subsystem."""
    print(f"\nCEED PROSPECTIVE FORECAST — issued {issue_year}")
    print("=" * 78)
    print(f"Initialised from {LAST_OBSERVED_YEAR} observations. Everything "
          f"below is out of sample.")
    print(f"Intervals are +/-{INTERVAL_SIGMA:.0f}x the HELD-OUT test RMSE, "
          f"not the in-sample fit.\n")

    for name in SYSTEMS:
        rows = [p for p in predictions if p.system == name]
        if not rows:
            continue
        head = rows[0]
        print(f"{name.upper()}  ({head.units})   "
              f"out-of-sample R2 = {head.test_r2:+.3f}  ->  {head.expectation}")
        print(f"  {'year':>6} {'predicted':>11} {'interval':>22} {'band':>16}")
        for p in rows:
            interval = f"[{p.lower:.2f}, {p.upper:.2f}]"
            print(f"  {p.year:>6} {p.value:>11.2f} {interval:>22} {p.band:>16}")
        print()


def trend_comparison(predictions: List[Prediction]) -> Dict[str, dict]:
    """Forecast trend against the observed 2010-2024 trend, per subsystem.

    Only meaningful for the monotonic subsystems. A linear trend fitted across
    part of an 11-year cycle is an artifact, so solar and magnetic are
    reported but flagged as not comparable.
    """
    cyclic = {'solar', 'magnetic'}
    out: Dict[str, dict] = {}
    for name in SYSTEMS:
        obs = ALL_OBSERVED[name]
        years = np.array(sorted(obs), dtype=float)
        values = np.array([obs[int(y)] for y in years], dtype=float)
        observed = float(np.polyfit(years, values, 1)[0])

        rows = [p for p in predictions if p.system == name]
        fy = np.array([p.year for p in rows], dtype=float)
        fv = np.array([p.value for p in rows], dtype=float)
        forecast = float(np.polyfit(fy, fv, 1)[0])

        out[name] = {
            "observed_per_year": observed,
            "forecast_per_year": forecast,
            "ratio": forecast / observed if observed else float('nan'),
            "comparable": name not in cyclic,
            "units": rows[0].units if rows else "",
        }
    return out


def falsification_criteria(predictions: List[Prediction]) -> List[str]:
    """The conditions under which this forecast counts as wrong.

    Written before resolution so they cannot be renegotiated afterwards.
    """
    by_system = {n: [p for p in predictions if p.system == n] for n in SYSTEMS}
    lines = []

    for name in SYSTEMS:
        rows = by_system[name]
        if not rows:
            continue
        expectation = rows[0].expectation
        if expectation == "expect skill":
            lines.append(
                f"{name}: FALSIFIED if fewer than half of the observed values "
                f"fall inside their intervals, or if a flat line at the 2024 "
                f"value scores a better RMSE over 2025-2030.")
        else:
            lines.append(
                f"{name}: this forecast is published EXPECTED TO FAIL "
                f"(out-of-sample R2 = {rows[0].test_r2:+.3f}). The prediction "
                f"under test is the failure itself. It would be a surprise — "
                f"and evidence against validation.py — if it verified well.")

    trends = trend_comparison(predictions)

    lines.append(
        "atmospheric, directional: the model has no ENSO term while a very "
        "strong El Nino is developing, so the 2026 and 2027 predictions are "
        "expected to run COLD. If they run WARM instead, the stated reason "
        "for the miss is wrong.")

    atm = trends['atmospheric']
    lines.append(
        f"atmospheric, quantitative: the forecast warms at "
        f"{atm['forecast_per_year']:+.3f} {atm['units']}/yr against an "
        f"observed 2010-2024 trend of {atm['observed_per_year']:+.3f}, i.e. "
        f"{atm['ratio']:.2f}x. Predicted miss is LOW. A high miss falsifies "
        f"the explanation even if the magnitude happens to fit.")

    oce = trends['oceanic']
    lines.append(
        f"oceanic, quantitative: the forecast accumulates at "
        f"{oce['forecast_per_year']:+.3f} {oce['units']}/yr against an "
        f"observed {oce['observed_per_year']:+.3f}, i.e. {oce['ratio']:.2f}x — "
        f"under-predicting by roughly {1/oce['ratio']:.1f}x. By 2030 that is "
        f"a shortfall of about "
        f"{(oce['observed_per_year']-oce['forecast_per_year'])*6:.1f} "
        f"{oce['units']}. This is the specific failure predicted; if oceanic "
        f"instead tracks observations, validation.py's -2.705 was misleading "
        f"and needs re-examining.")
    lines.append(
        "whole model: FALSIFIED as a forecasting tool if the mean out-of-"
        "sample R2 across subsystems over 2025-2030 is negative, i.e. the "
        "ensemble is beaten by persistence.")
    return lines


def to_record(predictions: List[Prediction], issue_year: int = 2026) -> dict:
    """Serialisable forecast record for archiving alongside the code."""
    return {
        "issued": issue_year,
        "initialised_from": LAST_OBSERVED_YEAR,
        "interval_sigma": INTERVAL_SIGMA,
        "interval_basis": "held-out test RMSE from simulation/validation.py",
        "model_has_enso": False,
        "predictions": [asdict(p) for p in predictions],
        "falsification": falsification_criteria(predictions),
        "resolution_sources": {
            "solar": "NOAA SWPC / NRCan Penticton F10.7 annual mean, sfu",
            "magnetic": "GFZ Potsdam Kp annual mean",
            "atmospheric": "NASA GISTEMP v4 + 0.19 K offset to 1850-1900 base",
            "oceanic": "NOAA/NCEI 0-700m OHC anomaly, 10^22 J",
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description="CEED prospective forecast (pre-registered)")
    parser.add_argument('--horizon', type=int, default=DEFAULT_HORIZON,
                        help=f'Years beyond {LAST_OBSERVED_YEAR}')
    parser.add_argument('--issue-year', type=int, default=2026)
    parser.add_argument('--json', type=str, default=None,
                        help='Write the machine-readable record here')
    args = parser.parse_args()

    print("Measuring out-of-sample skill before forecasting...")
    priors = skill_priors()
    for name in SYSTEMS:
        p = priors[name]
        print(f"  {name:<12} test R2 {p.test_r2:+.3f}  "
              f"test RMSE {p.test_rmse:.2f} {p.units:<10} -> {p.expectation}")

    predictions = run_forecast(horizon=args.horizon,
                               issue_year=args.issue_year,
                               priors=priors)
    print_forecast(predictions, issue_year=args.issue_year)

    print("-" * 78)
    print("FALSIFICATION CRITERIA — fixed at issue, not renegotiable")
    print("-" * 78)
    for line in falsification_criteria(predictions):
        print(f"  * {line}\n")

    if args.json:
        with open(args.json, 'w', encoding='utf-8') as fh:
            json.dump(to_record(predictions, args.issue_year), fh, indent=2)
        print(f"record written to {args.json}")


if __name__ == "__main__":
    main()
