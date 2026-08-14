# simulation/tipping.py

"""
Tipping elements: bistability, hysteresis, and commitment.

WHY THIS EXISTS
---------------
The convergence model has exactly one attractor. Integrate it from total
energies of 250, 500, 800 or 1500 and every trajectory lands on the same
state within 0.05 energy units. Its "state-dependent" parameters are pure
functions of (E_now, t) with no history argument, so a shock of any size
leaves no trace once it passes.

That makes it structurally unable to represent the thing the red team report
is about: a system pushed into a new state that does not come back. Two
observed cases motivated this module.

  Antarctic sea ice. The trend was slightly positive 1979-2015, then broke
  in September 2016. Record lows in 2023-2025 sit inside the new regime
  rather than being excursions from the old one, with "increased persistence
  in sea ice extent anomalies and a strongly reduced tendency to return to
  the mean state" (Comms Earth & Environment, 2025).

  Marine ice sheets. Tipping-element assessments find components that "can
  remain tipped even if the background climate falls back below the
  threshold" — the defining property of hysteresis.

WHAT THIS MODULE IS
-------------------
A dynamical primitive, not an ice sheet model. It is deliberately at the same
abstraction level as CEED_universal_model.py: the minimal system that exhibits
the three properties the convergence model lacks.

    dx/dt = (x - x^3 + F(t) + noise) / tau_fast     internal state
    dh/dt = (x - h) / tau_slow                       observable response

The cubic gives two stable branches separated by a fold; `tau_slow >>
tau_fast` separates the timescales.

The clearest physical reading is the shelf/sheet pair, which maps onto the
two variables exactly:

    x  fast   ice SHELF integrity / buttressing.  Floating, can disintegrate
              in weeks to months.  Larsen B lost ~3,250 km^2 in five weeks;
              the Thwaites eastern shelf is expected to go imminently.
    h  slow   grounded ice SHEET mass.  Responds over centuries to millennia.

They are coupled because the shelf buttresses the sheet: losing the shelf
removes the restraint, so a fast change in `x` commits `h` to a slow drawdown
that continues long after the trigger is gone. Shelf loss raises sea level by
almost nothing directly — it is already floating. It matters because of what
it stops holding back.

That is the whole reason for the timescale split. The event you can see is
fast and nearly harmless; the consequence it commits you to is slow and large.

THE POINT OF THE SPLIT
----------------------
`x` crosses the fold and commits long before `h` has visibly moved. The system
is decided while it still looks fine. `commitment_lag()` measures that window
directly, and it is the quantity a red team actually wants: not "what state am
I in" but "have I already chosen a state I cannot see yet".

Determinism: per the project rule, no random() inside any derivative. Noise is
pre-generated and interpolated, exactly as external events are in
convergence_model.py.
"""

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

# Fold (saddle-node bifurcation) points of dx/dt = x - x^3 + F.
# dx/dt = 0 and d/dx(x - x^3) = 0 give x = +/- 1/sqrt(3), F = -/+ 2/(3*sqrt(3)).
FOLD_X = 1.0 / np.sqrt(3.0)
FOLD_F = 2.0 / (3.0 * np.sqrt(3.0))  # ~0.3849


@dataclass
class TippingElement:
    """A fast-slow bistable element.

    Attributes:
        name: Human-readable label.
        tau_fast: Timescale of the internal state x [time].
        tau_slow: Timescale of the observable response h [time].
            Must exceed tau_fast for the commitment window to exist.
        x0: Initial internal state. Defaults to the lower stable branch.
        h0: Initial observable. Defaults to match x0.
    """
    name: str = "element"
    tau_fast: float = 1.0
    tau_slow: float = 50.0
    x0: float = -1.0
    h0: Optional[float] = None

    state: np.ndarray = field(init=False)

    def __post_init__(self):
        if self.tau_fast <= 0 or self.tau_slow <= 0:
            raise ValueError("timescales must be positive")
        if self.tau_slow <= self.tau_fast:
            raise ValueError(
                f"tau_slow ({self.tau_slow}) must exceed tau_fast "
                f"({self.tau_fast}); the commitment window is their separation")
        h0 = self.x0 if self.h0 is None else self.h0
        self.state = np.array([float(self.x0), float(h0)])

    # ── structure ────────────────────────────────────────────────────

    @staticmethod
    def fold_forcing() -> float:
        """|F| beyond which only one equilibrium remains."""
        return FOLD_F

    @staticmethod
    def equilibria(F: float) -> np.ndarray:
        """Real roots of x - x^3 + F = 0, sorted ascending.

        Returns 3 roots (two stable, one unstable) while |F| < FOLD_F,
        otherwise 1.
        """
        roots = np.roots([-1.0, 0.0, 1.0, F])
        real = np.sort(roots[np.abs(roots.imag) < 1e-9].real)
        return real

    @staticmethod
    def stable_equilibria(F: float) -> np.ndarray:
        """Equilibria where d/dx(x - x^3) = 1 - 3x^2 < 0."""
        eq = TippingElement.equilibria(F)
        return eq[(1.0 - 3.0 * eq ** 2) < 0]

    @staticmethod
    def is_bistable(F: float) -> bool:
        """True while two stable branches coexist."""
        return abs(F) < FOLD_F

    def is_committed(self, F: float, x: Optional[float] = None) -> bool:
        """Has the internal state crossed the fold toward the upper branch?

        This is the question the observable cannot answer. `h` may still sit
        near its old value while `x` has already passed the point of no
        return for the current forcing.
        """
        x = self.state[0] if x is None else x
        if F >= FOLD_F:
            return True          # upper branch is the only one left
        if F <= -FOLD_F:
            return False
        unstable = self.equilibria(F)
        if len(unstable) < 3:
            return x > 0
        return x > unstable[1]   # past the middle (unstable) root

    # ── dynamics ─────────────────────────────────────────────────────

    def derivative(self, state: Sequence[float], F: float,
                   noise: float = 0.0) -> np.ndarray:
        """RHS. Deterministic given (state, F, noise)."""
        x, h = state
        dx = (x - x ** 3 + F + noise) / self.tau_fast
        dh = (x - h) / self.tau_slow
        return np.array([dx, dh])

    def simulate(self,
                 forcing: Callable[[float], float],
                 years: float = 500.0,
                 dt: float = 0.1,
                 noise: Optional[np.ndarray] = None
                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Integrate with RK4.

        Args:
            forcing: F(t).
            years: Duration.
            dt: Timestep.
            noise: Pre-generated noise, one value per step. Generated OUTSIDE
                the derivative so the RHS stays deterministic.

        Returns:
            t, x, h — each of length steps + 1.
        """
        steps = int(round(years / dt))
        t = np.arange(steps + 1) * dt
        xs = np.empty(steps + 1)
        hs = np.empty(steps + 1)
        xs[0], hs[0] = self.state

        if noise is None:
            noise = np.zeros(steps + 1)
        elif len(noise) < steps + 1:
            raise ValueError(f"need >= {steps + 1} noise values, got {len(noise)}")

        s = self.state.copy()
        for i in range(1, steps + 1):
            ti = t[i - 1]
            n = noise[i - 1]
            k1 = self.derivative(s, forcing(ti), n)
            k2 = self.derivative(s + 0.5 * dt * k1, forcing(ti + 0.5 * dt), n)
            k3 = self.derivative(s + 0.5 * dt * k2, forcing(ti + 0.5 * dt), n)
            k4 = self.derivative(s + dt * k3, forcing(ti + dt), n)
            s = s + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            xs[i], hs[i] = s

        self.state = s
        return t, xs, hs

    # ── diagnostics ──────────────────────────────────────────────────

    def commitment_lag(self,
                       t: np.ndarray,
                       x: np.ndarray,
                       h: np.ndarray,
                       h_detect: float = 0.5) -> Optional[float]:
        """Time between committing and that becoming visible.

        Args:
            t, x, h: A trajectory from `simulate`.
            h_detect: Fraction of the full h excursion counted as "visible".

        Returns:
            Years between x crossing zero and h completing `h_detect` of its
            total change, or None if either never happens.
        """
        crossed = np.where(x > 0.0)[0]
        if len(crossed) == 0:
            return None
        t_commit = t[crossed[0]]

        span = h.max() - h.min()
        if span <= 0:
            return None
        threshold = h.min() + h_detect * span
        visible = np.where(h >= threshold)[0]
        if len(visible) == 0:
            return None
        return float(t[visible[0]] - t_commit)


# ── hysteresis ───────────────────────────────────────────────────────

def hysteresis_loop(element: Optional[TippingElement] = None,
                    F_max: float = 0.6,
                    ramp_years: float = 4000.0,
                    dt: float = 0.5
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Ramp forcing up past the fold, then symmetrically back down.

    The defining test: the return path does not retrace the outbound path.
    The system stays tipped well after the forcing that tipped it is gone.

    Returns:
        F_up, x_up, F_down, x_down
    """
    el = element or TippingElement(tau_fast=1.0, tau_slow=50.0, x0=-1.0)

    half = ramp_years / 2.0
    up = lambda t: -F_max + 2.0 * F_max * (t / half)
    t1, x1, _ = el.simulate(up, years=half, dt=dt)
    F_up = np.array([up(tt) for tt in t1])

    down = lambda t: F_max - 2.0 * F_max * (t / half)
    t2, x2, _ = el.simulate(down, years=half, dt=dt)
    F_down = np.array([down(tt) for tt in t2])

    return F_up, x1, F_down, x2


# ── early warning signals ────────────────────────────────────────────

def rolling_autocorrelation(series: np.ndarray, window: int) -> np.ndarray:
    """Lag-1 autocorrelation in a trailing window.

    Rises toward 1 as a system approaches a fold — recovery from
    perturbation slows down, so successive values resemble each other more.
    This is the standard early-warning signal, and the same statistic that
    exposed structure in the convergence model's hindcast residuals.
    """
    if window < 3:
        raise ValueError("window must be at least 3")
    out = np.full(len(series), np.nan)
    for i in range(window, len(series)):
        seg = series[i - window:i]
        seg = seg - seg.mean()
        denom = np.sum(seg ** 2)
        if denom > 0:
            out[i] = np.sum(seg[:-1] * seg[1:]) / denom
    return out


def rolling_variance(series: np.ndarray, window: int) -> np.ndarray:
    """Variance in a trailing window. Also rises approaching a fold."""
    if window < 3:
        raise ValueError("window must be at least 3")
    out = np.full(len(series), np.nan)
    for i in range(window, len(series)):
        out[i] = np.var(series[i - window:i])
    return out


def potential(x: float, F: float) -> float:
    """Potential U(x) whose gradient gives the fast dynamics.

    dx/dt = -dU/dx, so U(x) = -x^2/2 + x^4/4 - F*x.
    """
    return -x ** 2 / 2.0 + x ** 4 / 4.0 - F * x


def barrier_height(F: float) -> float:
    """Height of the potential barrier holding the system on the lower branch.

    This is what a shock has to clear to tip the system while the mean forcing
    is still below the fold. It collapses toward zero as F approaches the
    fold: 0.250 at F=0, 0.0252 at F=0.30, 0.00003 at F=0.384.

    A system described as "safely below threshold" needs a smaller and smaller
    kick the closer it gets.

    Returns:
        Barrier height, or 0.0 once no barrier exists (|F| >= FOLD_F).
    """
    if abs(F) >= FOLD_F:
        return 0.0
    eq = TippingElement.equilibria(F)
    if len(eq) < 3:
        return 0.0
    return potential(eq[1], F) - potential(eq[0], F)


@dataclass
class ForcingEvent:
    """A discrete forcing excursion — an atmospheric river, a storm, a pulse."""
    time: float
    magnitude: float
    duration: float = 5.0


def generate_events(rate: float,
                    years: float,
                    mean_magnitude: float,
                    rng: np.random.Generator,
                    heavy_tailed: bool = True,
                    duration: float = 5.0) -> List[ForcingEvent]:
    """Pre-generate a Poisson event train.

    Args:
        rate: Events per unit time.
        years: Horizon.
        mean_magnitude: Mean forcing excursion per event.
        rng: Seeded generator. Events are built BEFORE integration so the
            derivative stays deterministic.
        heavy_tailed: If True, magnitudes are exponential; if False, normal
            with standard deviation equal to the mean. Both then have the
            same mean AND the same variance, so any difference in outcome is
            tail shape alone.
        duration: How long each event holds.

    Returns:
        Events sorted by time.
    """
    n = rng.poisson(rate * years)
    times = np.sort(rng.uniform(0.0, years, n))
    if heavy_tailed:
        mags = rng.exponential(mean_magnitude, n)
    else:
        mags = rng.normal(mean_magnitude, mean_magnitude, n)
    return [ForcingEvent(float(t), float(m), duration)
            for t, m in zip(times, mags)]


def event_forcing(base: float,
                  events: Sequence[ForcingEvent]) -> Callable[[float], float]:
    """Build F(t) = base + the sum of any events active at t.

    Event arrays are hoisted out of the closure so each call is a vectorised
    comparison rather than a Python loop over every event. The forcing is
    evaluated three times per RK4 step, so the loop version dominated runtime.
    """
    if not events:
        return lambda t: base

    starts = np.array([e.time for e in events])
    ends = np.array([e.time + e.duration for e in events])
    mags = np.array([e.magnitude for e in events])

    def F(t: float) -> float:
        active = (starts <= t) & (t < ends)
        return base + float(mags[active].sum())
    return F


def escape_probability(base_forcing: float,
                       rate: float,
                       mean_magnitude: float,
                       trials: int = 200,
                       years: float = 500.0,
                       dt: float = 0.1,
                       heavy_tailed: bool = True,
                       seed: int = 0,
                       tau_fast: float = 1.0,
                       tau_slow: float = 50.0) -> float:
    """Fraction of trials that tip, with mean forcing held BELOW the fold.

    Bifurcation-induced tipping needs the forcing to cross the fold.
    Event-induced tipping does not: a single excursion can clear the barrier
    while the mean sits in the safe range. That is the atmospheric-river
    case — rare, brief, and sufficient.

    Raises:
        ValueError: if base_forcing is already at or beyond the fold, where
            the question is meaningless.
    """
    if abs(base_forcing) >= FOLD_F:
        raise ValueError(
            f"base_forcing {base_forcing} is already past the fold {FOLD_F}; "
            "no barrier remains to escape")

    rng = np.random.default_rng(seed)
    tipped = 0
    for _ in range(trials):
        events = generate_events(rate, years, mean_magnitude, rng,
                                 heavy_tailed=heavy_tailed)
        el = TippingElement(tau_fast=tau_fast, tau_slow=tau_slow, x0=-1.0)
        _, x, _ = el.simulate(event_forcing(base_forcing, events),
                              years=years, dt=dt)
        if np.any(x > 0.0):
            tipped += 1
    return tipped / trials


def recovery_rate(element: TippingElement, F: float,
                  branch: str = 'lower') -> float:
    """Linear recovery rate |d/dx(dx/dt)| at a stable equilibrium.

    Goes to zero at the fold — critical slowing down. The reciprocal is the
    recovery timescale, which diverges.

    Raises:
        ValueError: if the requested branch does not exist at this forcing.
    """
    stable = TippingElement.stable_equilibria(F)
    if len(stable) == 0:
        raise ValueError(f"no stable equilibrium at F={F}")
    if branch == 'lower':
        x_eq = stable[0]
    elif branch == 'upper':
        x_eq = stable[-1]
    else:
        raise ValueError("branch must be 'lower' or 'upper'")
    return abs((1.0 - 3.0 * x_eq ** 2) / element.tau_fast)


if __name__ == "__main__":
    if __package__ in (None, ""):
        # Allow `python simulation/tipping.py` as well as `-m simulation.tipping`
        import pathlib
        import sys
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

    print("CEED Tipping Element")
    print("=" * 66)

    el = TippingElement(name="demo", tau_fast=1.0, tau_slow=50.0)

    print(f"\nFold at |F| = {FOLD_F:.4f}, x = +/-{FOLD_X:.4f}")
    print("\nEquilibrium structure:")
    print(f"  {'F':>8} {'stable branches':>34} {'bistable':>10}")
    for F in (-0.5, -0.3, 0.0, 0.3, 0.38, 0.5):
        st = TippingElement.stable_equilibria(F)
        print(f"  {F:>8.2f} {np.array2string(st, precision=3):>34} "
              f"{str(TippingElement.is_bistable(F)):>10}")

    print("\nCritical slowing down (lower branch):")
    for F in (0.0, 0.2, 0.3, 0.37, 0.384):
        r = recovery_rate(el, F)
        print(f"  F={F:>6.3f}  recovery rate {r:>7.4f}  "
              f"recovery time {1.0 / r:>8.2f}")

    print("\n" + "-" * 66)
    print("COMMITMENT: the gap between deciding and showing")
    print("-" * 66)
    el = TippingElement(tau_fast=1.0, tau_slow=50.0, x0=-1.0)
    # Push just past the fold, then hold at a sub-threshold level forever.
    def forcing(t):
        return 0.45 if t < 100 else 0.0
    t, x, h = el.simulate(forcing, years=1200, dt=0.1)

    lag = el.commitment_lag(t, x, h)
    i_commit = np.where(x > 0)[0][0]
    print(f"  forcing removed at        t = 100")
    print(f"  internal state committed  t = {t[i_commit]:.1f}   "
          f"(h has moved {100*(h[i_commit]-h[0])/(h.max()-h[0]):.1f}% of its eventual change)")
    print(f"  response half-complete    t = {t[i_commit] + lag:.1f}")
    print(f"  commitment lag            = {lag:.1f} time units")
    print(f"\n  final x = {x[-1]:+.4f}   final h = {h[-1]:+.4f}")
    print("  forcing is long gone and the system has NOT returned.")

    print("\n" + "-" * 66)
    print("HYSTERESIS: the return path is not the outbound path")
    print("-" * 66)
    F_up, x_up, F_down, x_down = hysteresis_loop()
    up_switch = F_up[np.argmax(x_up > 0)]
    down_switch = F_down[np.argmax(x_down < 0)]
    print(f"  switches UP   at F = {up_switch:+.4f}")
    print(f"  switches DOWN at F = {down_switch:+.4f}")
    print(f"  loop width    = {up_switch - down_switch:.4f}  "
          f"(0 would mean no hysteresis)")

    print("\n" + "-" * 66)
    print("EVENT-INDUCED TIPPING: crossing without approaching")
    print("-" * 66)
    print("  Bifurcation tipping needs the mean forcing to reach the fold.")
    print("  A discrete excursion does not. Atmospheric rivers occupy ~3% of")
    print("  the time yet drive 40-80% of winter meltwater on peninsula")
    print("  shelves; rain on Thwaites reaches 30 mm in summer, 9 mm in")
    print("  winter. Hydrofracture then splits the shelf in weeks.\n")
    print(f"  {'base F':>8} {'barrier':>10} {'heavy tail':>12} {'thin tail':>11}")
    for F in (0.0, 0.15, 0.25):
        h = escape_probability(F, rate=0.05, mean_magnitude=0.10, trials=60,
                               years=400, dt=0.2, seed=5, heavy_tailed=True)
        th = escape_probability(F, rate=0.05, mean_magnitude=0.10, trials=60,
                                years=400, dt=0.2, seed=5, heavy_tailed=False)
        print(f"  {F:>8.2f} {barrier_height(F):>10.5f} {h:>11.0%} {th:>11.0%}")
    print("\n  At F=0 the barrier is at its maximum and the mean forcing is")
    print("  nowhere near the fold — and the system still tips. 'Below")
    print("  threshold' is not a safety claim once events are in the picture.")
    print("\n  Both event distributions have the SAME mean and variance; only")
    print("  the tail differs. The heavy tail dominates here because the fold")
    print("  sits 3.85x the typical event size away. Shrink that ratio and the")
    print("  ordering reverses — see E11 in legacy/README.md. Variance alone")
    print("  does not tell you the risk.")

    print("\n" + "-" * 66)
    print("CONTRAST: the convergence model under the same treatment")
    print("-" * 66)
    from simulation.convergence_model import (ConvergencePredictor,
                                              SystemParameters, SYSTEMS)
    p = SystemParameters()
    base = [p.E_initial[s] for s in SYSTEMS]
    finals = []
    for mult in (0.5, 1.0, 1.6, 3.0):
        pr = ConvergencePredictor(SystemParameters(
            E_initial=dict(zip(SYSTEMS, [mult * e for e in base]))))
        _, sol = pr.predict_convergence(years=200)
        finals.append(float(np.sum(sol[-1])))
        print(f"  start x{mult:<4} total {mult * sum(base):>7.1f}  "
              f"-> {finals[-1]:>7.1f}")
    print(f"  spread = {max(finals) - min(finals):.4f}  -> one attractor, no memory")
