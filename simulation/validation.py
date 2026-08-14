"""
Held-out validation for the CEED convergence model.

Everything the hindcast reported until now was measured on the same window
the parameters were tuned against. That is a fit report, not a test. This
module splits the record:

    train  2010-2019   parameters may see this
    test   2020-2024   parameters must never see this

and reports the gap between them. A model that fits the train window and
falls apart on the test window has memorised, not generalised.

Two properties make the split real:

1. ONE TRAJECTORY. The model is initialised from observations at 2010 and
   integrated straight through 2024. It is never re-initialised at 2020, so
   the test years are a free-running five-year forecast, not a warm start.

2. NO LEAKAGE INTO DERIVED CONSTANTS. The solar cycle amplitude is derived
   from the observed F10.7 range. Derived over the full record that range
   includes the 2024 maximum of 180 sfu, which sits in the test window — so
   the shipped constant leaks. `derive_solar_modulation` re-derives it from
   the train window alone.

Usage:
    python -m simulation.validation
    python -m simulation.validation --train-end 2018
"""

import argparse
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

from simulation.convergence_model import (
    SYSTEMS,
    AnthropogenicForcing,
    SystemParameters,
)
from simulation.hindcast import HindcastScore, run_hindcast
from simulation.unit_bridge import ALL_OBSERVED, SCALES

DEFAULT_START = 2010
DEFAULT_TRAIN_END = 2019
DEFAULT_END = 2024

# Systems whose fit the calibrated parameters actually control. A_0 and
# atm_fraction partition anthropogenic release between atmosphere and ocean;
# they have no meaningful lever on solar or magnetic.
CALIBRATION_TARGETS = ('atmospheric', 'oceanic')


# ── Derived constants, restricted to a window ────────────────────────

def derive_solar_modulation(params: SystemParameters,
                            first_year: int,
                            last_year: int) -> float:
    """Derive the solar source modulation from observations in a window.

    The solar subsystem is a first-order low-pass filter. Linearising about
    equilibrium gives an effective decay rate and hence an amplitude gain:

        lambda_eff = (lambda - alpha) + 2*gamma*E_eq
        gain       = 1/sqrt(1 + (omega/lambda_eff)^2)

    The source modulation must be pre-compensated for that gain to reproduce
    the observed half-amplitude:

        modulation = half_amplitude / (gain * E_eq)

    Args:
        params: Model parameters supplying alpha, lambda, gamma, E_eq.
        first_year: First observation year to use.
        last_year: Last observation year to use.

    Returns:
        Fractional source modulation.

    Raises:
        ValueError: if the window contains fewer than two observations.
    """
    obs = ALL_OBSERVED['solar']
    window = [v for yr, v in obs.items() if first_year <= yr <= last_year]
    if len(window) < 2:
        raise ValueError(
            f"need at least 2 solar observations in {first_year}-{last_year}, "
            f"got {len(window)}")

    E_eq = params.E_initial['solar']
    lam_eff = ((params.lambda_dissipation['solar']
                - params.alpha_retention['solar'])
               + 2.0 * params.gamma_nonlinear['solar'] * E_eq)
    omega = 2.0 * np.pi / 11.0
    gain = 1.0 / np.sqrt(1.0 + (omega / lam_eff) ** 2)

    # Observed half-amplitude, converted to energy units.
    half_amp_obs = (max(window) - min(window)) / 2.0
    half_amp_energy = half_amp_obs / SCALES['solar'].scale

    return half_amp_energy / (gain * E_eq)


# ── Calibration ──────────────────────────────────────────────────────

@dataclass
class CalibrationResult:
    """Outcome of fitting parameters on a training window."""
    params: SystemParameters
    A_0: float
    atm_fraction: float
    solar_modulation: float
    train_objective: float
    n_evaluations: int

    def summary(self) -> str:
        return (f"A_0={self.A_0:.3f}  atm_fraction={self.atm_fraction:.3f}  "
                f"solar_modulation={self.solar_modulation:.3f}  "
                f"train NRMSE={self.train_objective:.4f}  "
                f"({self.n_evaluations} evaluations)")


def _objective(params: SystemParameters,
               start_year: int,
               train_end: int) -> float:
    """Mean NRMSE over the calibration targets on the training window.

    NRMSE rather than RMSE so atmospheric (K) and oceanic (10^22 J) contribute
    comparably instead of the larger-numbered one dominating.
    """
    scores = run_hindcast(params,
                          systems=list(CALIBRATION_TARGETS),
                          start_year=start_year,
                          end_year=train_end)
    if not scores:
        return float('inf')
    values = [s.nrmse for s in scores.values()]
    if not all(np.isfinite(v) for v in values):
        return float('inf')
    return float(np.mean(values))


def calibrate(start_year: int = DEFAULT_START,
              train_end: int = DEFAULT_TRAIN_END,
              refine_rounds: int = 3,
              grid: int = 9) -> CalibrationResult:
    """Fit A_0 and atm_fraction on the training window only.

    Uses a deterministic coarse-to-fine grid search rather than a stochastic
    optimiser, so a given window always yields the same parameters and the
    result is reproducible without seeding.

    The solar modulation is DERIVED from the same window (not fitted), so no
    test-window observation reaches any parameter.

    Args:
        start_year: First year of the record (model initialised here).
        train_end: Last year the calibration may see.
        refine_rounds: Number of grid refinement passes.
        grid: Points per axis per pass.

    Returns:
        CalibrationResult holding the fitted parameters.
    """
    base = SystemParameters()
    solar_mod = derive_solar_modulation(base, start_year, train_end)

    def build(A_0: float, atm_fraction: float) -> SystemParameters:
        return SystemParameters(
            solar_modulation=solar_mod,
            anthropogenic=AnthropogenicForcing(
                A_0=A_0, atm_fraction=atm_fraction),
        )

    A_lo, A_hi = 0.5, 6.0
    f_lo, f_hi = 0.05, 0.60

    best = (float('inf'), None, None)
    evaluations = 0

    for _ in range(refine_rounds):
        for A_0 in np.linspace(A_lo, A_hi, grid):
            for frac in np.linspace(f_lo, f_hi, grid):
                score = _objective(build(A_0, frac), start_year, train_end)
                evaluations += 1
                if score < best[0]:
                    best = (score, float(A_0), float(frac))

        # Shrink the bracket around the incumbent for the next pass.
        _, A_best, f_best = best
        A_span = (A_hi - A_lo) / (grid - 1)
        f_span = (f_hi - f_lo) / (grid - 1)
        A_lo, A_hi = A_best - A_span, A_best + A_span
        f_lo, f_hi = max(0.01, f_best - f_span), min(0.99, f_best + f_span)

    objective, A_best, f_best = best
    return CalibrationResult(
        params=build(A_best, f_best),
        A_0=A_best,
        atm_fraction=f_best,
        solar_modulation=solar_mod,
        train_objective=objective,
        n_evaluations=evaluations,
    )


# ── Split evaluation ─────────────────────────────────────────────────

@dataclass
class SplitResult:
    """Train and test scores for one parameter set."""
    label: str
    train: Dict[str, HindcastScore]
    test: Dict[str, HindcastScore]
    train_years: Tuple[int, int]
    test_years: Tuple[int, int]

    def mean_r2(self, which: str = 'test') -> float:
        scores = self.test if which == 'test' else self.train
        if not scores:
            return float('nan')
        return float(np.mean([s.r_squared for s in scores.values()]))

    def gap(self, system: str) -> float:
        """Train R² minus test R². Large positive = overfit."""
        return self.train[system].r_squared - self.test[system].r_squared


def evaluate_split(params: SystemParameters,
                   label: str,
                   start_year: int = DEFAULT_START,
                   train_end: int = DEFAULT_TRAIN_END,
                   end_year: int = DEFAULT_END,
                   systems: Optional[List[str]] = None) -> SplitResult:
    """Score one parameter set on the train and test windows.

    Both windows come from a SINGLE trajectory initialised at start_year, so
    the test window is a free-running forecast.
    """
    train = run_hindcast(params, systems=systems,
                         start_year=start_year, end_year=end_year,
                         score_from=start_year + 1, score_to=train_end)
    test = run_hindcast(params, systems=systems,
                        start_year=start_year, end_year=end_year,
                        score_from=train_end + 1, score_to=end_year)
    return SplitResult(
        label=label,
        train=train, test=test,
        train_years=(start_year + 1, train_end),
        test_years=(train_end + 1, end_year),
    )


def print_split(result: SplitResult):
    """Print a train/test comparison table."""
    tr0, tr1 = result.train_years
    te0, te1 = result.test_years

    print(f"\n{result.label}")
    print("=" * 78)
    print(f"{'System':<14s} {'train R²':>9s} {'test R²':>9s} {'gap':>8s} "
          f"{'test RMSE':>12s} {'test grade':>11s}")
    print(f"{'':<14s} {f'{tr0}-{tr1}':>9s} {f'{te0}-{te1}':>9s}")
    print("-" * 78)

    for name in result.test:
        tr = result.train[name]
        te = result.test[name]
        print(f"{name:<14s} {tr.r_squared:>+9.3f} {te.r_squared:>+9.3f} "
              f"{result.gap(name):>+8.3f} "
              f"{te.rmse:>8.2f} {te.units:<3s} {te.grade():>11s}")

    print("-" * 78)
    print(f"{'Mean':<14s} {result.mean_r2('train'):>+9.3f} "
          f"{result.mean_r2('test'):>+9.3f} "
          f"{result.mean_r2('train') - result.mean_r2('test'):>+8.3f}")


def run_validation(start_year: int = DEFAULT_START,
                   train_end: int = DEFAULT_TRAIN_END,
                   end_year: int = DEFAULT_END
                   ) -> Tuple[CalibrationResult, SplitResult, SplitResult]:
    """Calibrate on train, then score shipped and refitted parameters.

    Returns:
        (calibration, shipped_split, refitted_split)
    """
    print(f"Calibrating on {start_year}-{train_end}, "
          f"holding out {train_end + 1}-{end_year}...")
    calibration = calibrate(start_year=start_year, train_end=train_end)
    print(f"  {calibration.summary()}")

    shipped = SystemParameters(anthropogenic=AnthropogenicForcing())
    shipped_split = evaluate_split(
        shipped, "SHIPPED PARAMETERS (tuned on the full record — leaks)",
        start_year, train_end, end_year)

    refitted_split = evaluate_split(
        calibration.params,
        f"REFITTED ON {start_year}-{train_end} ONLY (no leakage)",
        start_year, train_end, end_year)

    return calibration, shipped_split, refitted_split


def main():
    parser = argparse.ArgumentParser(
        description='CEED held-out validation (train/test split)')
    parser.add_argument('--start', type=int, default=DEFAULT_START,
                        help=f'First year (default {DEFAULT_START})')
    parser.add_argument('--train-end', type=int, default=DEFAULT_TRAIN_END,
                        help=f'Last training year (default {DEFAULT_TRAIN_END})')
    parser.add_argument('--end', type=int, default=DEFAULT_END,
                        help=f'Last year (default {DEFAULT_END})')
    args = parser.parse_args()

    calibration, shipped, refitted = run_validation(
        args.start, args.train_end, args.end)

    print_split(shipped)
    print_split(refitted)

    print("\n" + "=" * 78)
    print("READING THIS")
    print("-" * 78)
    print("  gap = train R² - test R².  Large positive means the fit does not")
    print("  carry forward.  Negative means the test window happened to be")
    print("  easier, which is not evidence of skill either.")
    print()
    print(f"  The test window is {args.end - args.train_end} years and the solar")
    print("  cycle is 11, so this split cannot test a full cycle. Treat it as")
    print("  a check on drift and bias, not a verdict on cycle physics.")
    print()
    print("  A negative test R² means the model is beaten by a flat line drawn")
    print("  through the test window's own mean.")


if __name__ == "__main__":
    main()
