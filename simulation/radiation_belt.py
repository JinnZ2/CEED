"""Radiation-belt barrier: an L slot for the magnetic subsystem, crossed by events.

Co-authored, human + Claude, 2026-10-28. Motivated by the May 2024 superstorm:
Arase saw multi-MeV electrons inside the 2014 "impenetrable barrier" (Miyoshi et
al. 2026, GRL, doi:10.1029/2026GL125879) and CALET saw a new population at
L 2.2-3.2 that lasted months (arXiv:2602.03990). See Docs/arase-may2024/.

Two proposals wired in here:

1. The magnetic subsystem has no slot for WHERE energy lands. Kp is a global
   index. This module adds two stores, outside the barrier (outer belt) and
   inside it, driven by the model's own magnetic trajectory converted to Kp.

2. Event-induced crossing (the shape of tipping.escape_probability). The
   barrier is not a wall: it is the largest storm in the window it was
   characterized on (`characterized_barrier`). The mean Kp can sit far below
   it while one storm from the tail crosses. Storms are pre-generated with
   tipping.generate_events, so nothing random happens inside the derivative.

Units: stores are in ARBITRARY units (ILLUSTRATIVE). Driver is Kp. Time in days
internally; the public entry takes model years and converts.

Kp ceiling: Kp is defined on 0..9 (KP_MAX), so the driver saturates. Two storms
far apart in size both read 9. A barrier characterized on a window that already
reached Kp 9 cannot be crossed by this driver at all: the instrument has no
reading above its own ceiling. `BeltResult.saturated_days` counts it, and
`characterized_barrier` returns AT_CEILING rather than a number in that case.
An unsaturated storm index (Dst, SYM-H) is not in CEED; that is the gap.

Constants:
  INNER_EFOLD_DAYS = 41.0  CALET e-folding for >1.5 MeV at L 2.8-3.2, as read
                           from the arXiv abstract page. VERIFY against the paper.
  Everything else (injection rate, outer-belt fill and decay, storm rate and
  size) is ILLUSTRATIVE: no literature value is carried here, and none should
  be read out of this module as a prediction. Not calibrated.

Not wired into the ODE: belt dynamics run on days, the convergence model on
months. The belt reads the magnetic trajectory; it does not feed back.
"""
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

from simulation.tipping import ForcingEvent, event_forcing, generate_events

DAYS_PER_YEAR = 365.25
KP_MAX = 9.0                    # top of the Kp scale by definition

# CALET (arXiv:2602.03990), >1.5 MeV, L 2.8-3.2. VERIFY.
INNER_EFOLD_DAYS = 41.0

# ILLUSTRATIVE: not traceable to literature. Do not quote as physics.
OUTER_FILL_PER_KP_DAY = 1.0     # outer-belt source per Kp per day
OUTER_EFOLD_DAYS = 10.0         # outer-belt loss timescale
INJECT_PER_KP_DAY = 0.2         # fraction of outer store moved inward per Kp of excess per day


@dataclass
class BeltResult:
    """Output of a belt run. `inner` is the store inside the barrier."""
    t_days: np.ndarray
    kp: np.ndarray
    outer: np.ndarray
    inner: np.ndarray
    barrier_kp: float
    crossing_days: List[float]
    saturated_days: int

    @property
    def crossed(self) -> bool:
        return bool(self.crossing_days)


def characterized_barrier(window_kp: Sequence[float]) -> float:
    """The empirical barrier: the largest driver seen in the window it was mapped on.

    A boundary characterized on a window is a reading of that window. A driver
    larger than anything in the window is outside the conditions the boundary
    was measured under (error type E2, stale conditions, if read as a wall).
    """
    w = np.clip(np.asarray(window_kp, dtype=float), 0.0, KP_MAX)
    if w.size == 0:
        raise ValueError("cannot characterize a barrier on an empty window")
    return float(np.max(w))


def barrier_reading(barrier_kp: float) -> str:
    """AT_CEILING if the barrier sits at Kp 9: no Kp reading can exceed it."""
    return "AT_CEILING" if barrier_kp >= KP_MAX else "BELOW_CEILING"


def storm_train(rate_per_year: float, years: float, mean_kp: float,
                seed: int = 0, heavy_tailed: bool = True,
                duration_days: float = 2.0) -> List[ForcingEvent]:
    """Pre-generate storms (Kp excursions) with tipping.generate_events.

    Times and durations are returned in DAYS. Same mean and variance in both
    tail shapes, as in tipping.py, so a difference in crossings is tail shape.
    """
    rng = np.random.default_rng(seed)
    ev = generate_events(rate_per_year, years, mean_kp, rng,
                         heavy_tailed=heavy_tailed, duration=1.0)
    return [ForcingEvent(e.time * DAYS_PER_YEAR, max(e.magnitude, 0.0), duration_days)
            for e in ev]


def run_belt(base_kp: np.ndarray, t_days: np.ndarray,
             storms: Sequence[ForcingEvent], barrier_kp: float,
             inner_efold_days: float = INNER_EFOLD_DAYS,
             outer0: float = 0.0, inner0: float = 0.0) -> BeltResult:
    """Integrate the two stores on a daily grid (explicit Euler, dt = grid step).

    d outer/dt = FILL*kp - outer/OUTER_EFOLD - inj
    d inner/dt = inj - inner/inner_efold
    inj        = INJECT * max(kp - barrier, 0) * outer   (zero below the barrier)

    kp(t) = base_kp(t) + sum of active storms. Deterministic: storms are inputs.
    """
    t_days = np.asarray(t_days, dtype=float)
    base_kp = np.asarray(base_kp, dtype=float)
    if t_days.shape != base_kp.shape or t_days.size < 2:
        raise ValueError("t_days and base_kp must be equal-length arrays, length >= 2")
    F = event_forcing(0.0, list(storms))
    raw = base_kp + np.array([F(t) for t in t_days])
    kp = np.clip(raw, 0.0, KP_MAX)
    saturated = int(np.sum(raw > KP_MAX))
    outer = np.zeros_like(t_days)
    inner = np.zeros_like(t_days)
    outer[0], inner[0] = outer0, inner0
    crossings: List[float] = []
    above_prev = False
    for k in range(1, t_days.size):
        dt = t_days[k] - t_days[k - 1]
        excess = max(kp[k - 1] - barrier_kp, 0.0)
        inj = INJECT_PER_KP_DAY * excess * outer[k - 1]
        outer[k] = outer[k - 1] + dt * (OUTER_FILL_PER_KP_DAY * kp[k - 1]
                                        - outer[k - 1] / OUTER_EFOLD_DAYS - inj)
        inner[k] = inner[k - 1] + dt * (inj - inner[k - 1] / inner_efold_days)
        above = excess > 0.0
        if above and not above_prev:
            crossings.append(float(t_days[k - 1]))
        above_prev = above
    return BeltResult(t_days, kp, outer, inner, barrier_kp, crossings, saturated)


def belt_from_magnetic(t_years: np.ndarray, magnetic_E: np.ndarray,
                       storms: Sequence[ForcingEvent], barrier_kp: float,
                       dt_days: float = 1.0) -> BeltResult:
    """Drive the belt from the convergence model's magnetic column.

    Converts magnetic energy to Kp with unit_bridge, interpolates onto a daily
    grid, and runs `run_belt`. This is the wiring point; the ODE is unchanged.
    """
    from simulation.unit_bridge import energy_to_observations
    t_years = np.asarray(t_years, dtype=float)
    kp_model = energy_to_observations('magnetic', np.asarray(magnetic_E, dtype=float))
    t_days = np.arange(0.0, t_years[-1] * DAYS_PER_YEAR + dt_days / 2, dt_days)
    base = np.interp(t_days, t_years * DAYS_PER_YEAR, kp_model)
    return run_belt(base, t_days, storms, barrier_kp)


def main() -> None:
    years = 10.0
    t_days = np.arange(0.0, years * DAYS_PER_YEAR, 1.0)
    base = np.full_like(t_days, 2.0)            # quiet mean Kp, ILLUSTRATIVE
    window = storm_train(4.0, 2.0, 1.5, seed=1)  # storms the barrier was mapped on
    barrier = characterized_barrier([2.0 + s.magnitude for s in window])
    print(f"barrier (largest Kp in the 2-year characterization window): "
          f"{barrier:.2f} [{barrier_reading(barrier)}]")
    print(f"mean Kp {base.mean():.2f}: far below the barrier\n")
    print(f"{'tail':<6} {'runs crossing':>14} {'mean crossings':>15} {'runs reaching Kp 9':>19}")
    n = 200
    for heavy in (True, False):
        hit = cnt = sat = 0
        for seed in range(n):
            storms = storm_train(4.0, years, 1.5, seed=100 + seed, heavy_tailed=heavy)
            r = run_belt(base, t_days, storms, barrier)
            hit += r.crossed
            cnt += len(r.crossing_days)
            sat += r.saturated_days > 0
        print(f"{'heavy' if heavy else 'thin':<6} {hit / n:>14.3f} {cnt / n:>15.2f} {sat / n:>19.3f}")
    print(f"\n{n} ten-year runs per tail; same mean and variance per storm.")
    print("\nInner-store units are ILLUSTRATIVE. Only the inner e-folding "
          f"({INNER_EFOLD_DAYS:.0f} d, CALET, VERIFY) is carried from literature.")


if __name__ == "__main__":
    main()
