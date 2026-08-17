# simulation/cascade.py

"""
Tipping cascades: coupled bistable elements, where one tipping shifts another's
threshold.

WHY THIS EXISTS
---------------
CEED stands for Cascading Energetic Event Disruption. Until now nothing in the
repository modelled a cascade. The convergence model has one attractor, so
nothing in it can tip at all; `tipping.py` supplies a single element that can.
This module couples them, which is the actual thesis.

THE STRUCTURE
-------------
For N elements, each with the fast-slow form from `tipping.py`:

    dx_i/dt = (x_i - x_i^3 + F_i(t) + sum_j C_ij * phi_j + noise) / tau_fast_i
    dh_i/dt = (x_i - h_i) / tau_slow_i
    phi_j   = (h_j + 1) / 2                       element j's tipped fraction

`C_ij` is the effect of element j on element i's effective forcing.
Positive is destabilising (j tipping pushes i toward its fold), negative is
stabilising. The matrix is NOT assumed symmetric — Thwaites influencing Ross
is not the same relationship as Ross influencing Thwaites.

Coupling runs through `phi` rather than raw `h` so an untouched network exerts
no influence at all. Coupling to `h` directly would give every element a
spurious shove at t=0, and would invert stabilising links, since a negative
coefficient times h = -1 is a positive push.

COUPLING RUNS THROUGH h, NOT x
------------------------------
This is the important design choice. Element j affects element i through j's
SLOW variable, so the influence only arrives as j actually responds. A shelf
disintegrating does not instantly alter a neighbouring system; the grounded
ice drawdown it commits to does, over centuries.

That makes cascades rate- and delay-dependent, which the literature finds is
decisive: interacting-tipping-element assessments report that whether one
element destabilises or stabilises another is "contingent upon the rates of
and delays between their collapse". Coupling through x would erase exactly the
effect that matters.

WHAT THIS IS NOT
----------------
Not a calibrated ice sheet model. The Thwaites/Ross configuration in
`__main__` is an illustrative parameterisation chosen to expose the asymmetry
between proximity-to-threshold and size-of-consequence. Do not read numbers
out of it as predictions.
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    from simulation.tipping import FOLD_F
except ModuleNotFoundError:  # running as `python simulation/cascade.py`
    import pathlib
    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
    from simulation.tipping import FOLD_F


@dataclass
class NetworkElement:
    """One node: a tipping element plus what tipping it would cost."""
    name: str
    tau_fast: float = 1.0
    tau_slow: float = 50.0
    x0: float = -1.0
    consequence: float = 1.0   # arbitrary units of "what is lost if it tips"

    def __post_init__(self):
        if self.tau_slow <= self.tau_fast:
            raise ValueError(
                f"{self.name}: tau_slow must exceed tau_fast")


@dataclass
class TippingNetwork:
    """N coupled fast-slow bistable elements.

    Attributes:
        elements: Nodes, in the order the coupling matrix indexes them.
        coupling: (N, N) matrix. coupling[i, j] is the effect of element j's
            slow state on element i's forcing. Diagonal must be zero — an
            element's self-feedback already lives in its cubic.
    """
    elements: List[NetworkElement]
    coupling: np.ndarray

    state: np.ndarray = field(init=False)

    def __post_init__(self):
        n = len(self.elements)
        self.coupling = np.asarray(self.coupling, dtype=float)
        if self.coupling.shape != (n, n):
            raise ValueError(
                f"coupling must be ({n},{n}), got {self.coupling.shape}")
        if np.any(np.diag(self.coupling) != 0.0):
            raise ValueError("coupling diagonal must be zero")
        x = np.array([e.x0 for e in self.elements], dtype=float)
        self.state = np.concatenate([x, x.copy()])

    @property
    def n(self) -> int:
        return len(self.elements)

    @property
    def names(self) -> List[str]:
        return [e.name for e in self.elements]

    def split(self, state: Optional[np.ndarray] = None
              ) -> Tuple[np.ndarray, np.ndarray]:
        """Return (x, h) views of a flat state vector."""
        s = self.state if state is None else np.asarray(state, dtype=float)
        return s[:self.n], s[self.n:]

    # ── dynamics ─────────────────────────────────────────────────────

    def derivative(self, state: Sequence[float], forcings: np.ndarray,
                   noise: Optional[np.ndarray] = None) -> np.ndarray:
        """RHS. Deterministic given (state, forcings, noise)."""
        x, h = self.split(np.asarray(state, dtype=float))
        if noise is None:
            noise = np.zeros(self.n)

        # Influence arrives through the slow variables, expressed as each
        # element's TIPPED FRACTION rather than its raw h.
        #
        # h runs from -1 (intact) to +1 (gone), so coupling to h directly
        # would make an untouched network exert influence at t=0 — and would
        # invert the sign of a stabilising link, since a negative coefficient
        # times h = -1 is a positive shove. phi = (h + 1) / 2 is 0 when
        # intact and 1 when fully tipped, so influence is zero until
        # something actually tips and a negative coefficient is protective.
        phi = (h + 1.0) / 2.0
        influence = self.coupling @ phi

        tau_f = np.array([e.tau_fast for e in self.elements])
        tau_s = np.array([e.tau_slow for e in self.elements])

        dx = (x - x ** 3 + forcings + influence + noise) / tau_f
        dh = (x - h) / tau_s
        return np.concatenate([dx, dh])

    def simulate(self,
                 forcing: Callable[[float], Sequence[float]],
                 years: float = 2000.0,
                 dt: float = 0.5,
                 noise: Optional[np.ndarray] = None
                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Integrate with RK4.

        Args:
            forcing: F(t) -> array of N forcings.
            years: Duration.
            dt: Timestep.
            noise: Pre-generated, shape (steps + 1, N). Never drawn inside
                the derivative, per the project rule.

        Returns:
            t, x, h — x and h have shape (steps + 1, N).
        """
        steps = int(round(years / dt))
        t = np.arange(steps + 1) * dt
        xs = np.empty((steps + 1, self.n))
        hs = np.empty((steps + 1, self.n))
        xs[0], hs[0] = self.split()

        if noise is None:
            noise = np.zeros((steps + 1, self.n))
        elif noise.shape[0] < steps + 1:
            raise ValueError(
                f"need >= {steps + 1} noise rows, got {noise.shape[0]}")

        s = self.state.copy()
        for i in range(1, steps + 1):
            ti, nz = t[i - 1], noise[i - 1]
            F1 = np.asarray(forcing(ti), dtype=float)
            F2 = np.asarray(forcing(ti + 0.5 * dt), dtype=float)
            F3 = np.asarray(forcing(ti + dt), dtype=float)
            k1 = self.derivative(s, F1, nz)
            k2 = self.derivative(s + 0.5 * dt * k1, F2, nz)
            k3 = self.derivative(s + 0.5 * dt * k2, F2, nz)
            k4 = self.derivative(s + dt * k3, F3, nz)
            s = s + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

            if not np.all(np.isfinite(s)):
                # The cubic diverges under RK4 once dt approaches the fastest
                # timescale. Fail loudly rather than propagate NaN into the
                # diagnostics, where "never tipped" and "blew up" look alike.
                raise FloatingPointError(
                    f"integration diverged at t={t[i]:.3f} with dt={dt}. "
                    f"Reduce dt below min(tau_fast)="
                    f"{min(e.tau_fast for e in self.elements)}.")

            xs[i], hs[i] = self.split(s)

        self.state = s
        return t, xs, hs

    # ── diagnostics ──────────────────────────────────────────────────

    def tipped(self, x_row: np.ndarray) -> np.ndarray:
        """Boolean mask of which elements are on the upper branch."""
        return np.asarray(x_row) > 0.0

    def tip_times(self, t: np.ndarray, x: np.ndarray
                  ) -> Dict[str, Optional[float]]:
        """First time each element crosses to the upper branch."""
        out: Dict[str, Optional[float]] = {}
        for i, name in enumerate(self.names):
            crossed = np.where(x[:, i] > 0.0)[0]
            out[name] = float(t[crossed[0]]) if len(crossed) else None
        return out

    def consequence_realised(self, h_row: np.ndarray) -> float:
        """Sum of consequence weighted by how far each slow variable has moved.

        h runs from -1 (intact) to +1 (fully gone), so (h + 1) / 2 is the
        fraction realised. Clipped to [0, 1]: destabilising coupling shifts an
        element's equilibrium beyond +1, which means "fully tipped", not
        "more than fully tipped".
        """
        weights = np.array([e.consequence for e in self.elements])
        fraction = np.clip((np.asarray(h_row) + 1.0) / 2.0, 0.0, 1.0)
        return float(np.sum(weights * fraction))

    def basin_margin(self, x_row: Optional[np.ndarray] = None,
                     forcings: Optional[np.ndarray] = None) -> np.ndarray:
        """Distance from each element's state to its basin boundary.

        The boundary is the unstable (middle) equilibrium. Small margin means
        near-tipping. This is the meaningful "how close is it" measure —
        forcing margin alone cannot distinguish two elements sitting at
        different points within the same basin.

        Returns:
            Signed distance per element. Negative means already past.
        """
        x = self.split()[0] if x_row is None else np.asarray(x_row, dtype=float)
        F = np.zeros(self.n) if forcings is None else np.asarray(forcings)

        margins = np.empty(self.n)
        for i in range(self.n):
            roots = np.roots([-1.0, 0.0, 1.0, F[i]])
            real = np.sort(roots[np.abs(roots.imag) < 1e-9].real)
            boundary = real[1] if len(real) == 3 else (
                -np.inf if F[i] >= FOLD_F else np.inf)
            margins[i] = boundary - x[i]
        return margins

    def cascade_from(self,
                     trigger: str,
                     base_forcing: Optional[np.ndarray] = None,
                     push: float = 0.45,
                     push_years: float = 100.0,
                     years: float = 4000.0,
                     dt: float = 0.5
                     ) -> Dict[str, Optional[float]]:
        """Force one element past its fold, release, and see what follows.

        The trigger is pushed for `push_years` then the forcing returns to
        `base_forcing`. Anything that tips afterwards did so because of the
        coupling, not because it was being forced.

        Args:
            trigger: Name of the element to push.
            base_forcing: Background forcing on every element. Default zeros.
            push: Forcing applied to the trigger. Must exceed the fold.
            push_years: How long the push lasts.
            years: Total integration.
            dt: Timestep.

        Returns:
            Tip time per element, or None for elements that never tipped.

        Raises:
            ValueError: if `trigger` is not in the network.
        """
        if trigger not in self.names:
            raise ValueError(f"unknown element {trigger!r}; have {self.names}")
        idx = self.names.index(trigger)
        base = (np.zeros(self.n) if base_forcing is None
                else np.asarray(base_forcing, dtype=float))

        def forcing(t):
            F = base.copy()
            if t < push_years:
                F[idx] = push
            return F

        # Reset to initial conditions before running.
        x0 = np.array([e.x0 for e in self.elements], dtype=float)
        self.state = np.concatenate([x0, x0.copy()])

        t, x, h = self.simulate(forcing, years=years, dt=dt)
        return self.tip_times(t, x)

    def domino_threshold(self,
                         trigger: str,
                         target: str,
                         lo: float = 0.0,
                         hi: float = 2.0,
                         tolerance: float = 1e-3,
                         **kwargs) -> Optional[float]:
        """Smallest coupling strength at which tripping `trigger` tips `target`.

        Bisects on coupling[target, trigger], restoring it afterwards.

        Returns:
            The threshold, or None if `target` never tips even at `hi`.
        """
        i, j = self.names.index(target), self.names.index(trigger)
        original = self.coupling[i, j]

        def tips(strength: float) -> bool:
            self.coupling[i, j] = strength
            return self.cascade_from(trigger, **kwargs)[target] is not None

        try:
            if not tips(hi):
                return None
            if tips(lo):
                return lo
            while hi - lo > tolerance:
                mid = 0.5 * (lo + hi)
                if tips(mid):
                    hi = mid
                else:
                    lo = mid
            return 0.5 * (lo + hi)
        finally:
            self.coupling[i, j] = original


# ── an illustrative two-element configuration ────────────────────────

def thwaites_ross_network(coupling_strength: float = 0.55,
                          back_coupling: float = 0.0) -> TippingNetwork:
    """A two-element network showing the proximity/consequence asymmetry.

    ILLUSTRATIVE ONLY — not calibrated, not a prediction. The parameters
    encode a qualitative structure reported in the literature:

      Thwaites   warm cavity, modified Circumpolar Deep Water warming since
                 2020, "largest changes of any ice-ocean system in
                 Antarctica". Near its fold. Fast. Modest direct consequence.
      Ross       cold cavity, currently stable, far from its fold. Very slow.
                 Buttresses a catchment holding ~11.6 m of sea level
                 equivalent, so the consequence is an order of magnitude
                 larger.

    The element closest to tipping is not the one that matters.

    Args:
        coupling_strength: Thwaites -> Ross influence. Positive destabilises.
        back_coupling: Ross -> Thwaites influence.
    """
    elements = [
        NetworkElement(name="thwaites", tau_fast=1.0, tau_slow=40.0,
                       x0=-0.6, consequence=0.65),
        NetworkElement(name="ross", tau_fast=5.0, tau_slow=400.0,
                       x0=-1.0, consequence=11.6),
    ]
    coupling = np.array([
        [0.0,               back_coupling],   # onto thwaites
        [coupling_strength, 0.0],             # onto ross
    ])
    return TippingNetwork(elements=elements, coupling=coupling)


if __name__ == "__main__":
    print("CEED Tipping Cascade")
    print("=" * 70)
    print("Illustrative configuration — not calibrated, not a prediction.\n")

    net = thwaites_ross_network()
    margins = net.basin_margin()
    print(f"{'element':<10} {'x0':>6} {'tau_fast':>9} {'tau_slow':>9} "
          f"{'consequence':>12} {'basin margin':>13}")
    for e, m in zip(net.elements, margins):
        print(f"{e.name:<10} {e.x0:>6.2f} {e.tau_fast:>9.1f} "
              f"{e.tau_slow:>9.1f} {e.consequence:>12.2f} {m:>13.3f}")
    print("  basin margin = distance to the unstable root; smaller is closer "
          "to tipping")

    print("\n" + "-" * 70)
    print("1. UNCOUPLED: tip Thwaites, nothing else should follow")
    print("-" * 70)
    solo = thwaites_ross_network(coupling_strength=0.0)
    print("  tip times:", solo.cascade_from("thwaites"))

    print("\n" + "-" * 70)
    print("2. COUPLED: same push, coupling on")
    print("-" * 70)
    net = thwaites_ross_network(coupling_strength=0.55)
    times = net.cascade_from("thwaites")
    for name, tt in times.items():
        print(f"  {name:<10} tips at t = {tt if tt is None else round(tt, 1)}")
    if times["ross"] is not None:
        print(f"\n  Ross tips {times['ross'] - times['thwaites']:.0f} time units "
              f"after Thwaites, with no forcing of its own.")

    print("\n" + "-" * 70)
    print("3. DOMINO THRESHOLD: how much coupling does it take?")
    print("-" * 70)
    probe = thwaites_ross_network(coupling_strength=0.0)
    thr = probe.domino_threshold("thwaites", "ross")
    print(f"  minimum Thwaites->Ross coupling for a cascade: {thr:.4f}"
          if thr else "  no cascade at any tested strength")

    print("\n" + "-" * 70)
    print("4. SIGN AND TIMING BOTH MATTER")
    print("-" * 70)
    f = lambda v: "never" if v is None else f"{v:.0f}"

    print("  (a) forward link Thwaites->Ross: sign decides whether it spreads")
    print(f"  {'coupling':>16} {'thwaites':>12} {'ross':>12}")
    for c in (0.55, 0.30, -0.55):
        n = thwaites_ross_network(coupling_strength=c)
        tt = n.cascade_from("thwaites")
        print(f"  {c:>16.2f} {f(tt['thwaites']):>12} {f(tt['ross']):>12}")

    print("\n  (b) back link Ross->Thwaites: arrives too late to protect it")
    print(f"  {'back coupling':>16} {'thwaites':>12} {'ross':>12}")
    for back in (0.0, -0.6, -2.0):
        n = thwaites_ross_network(coupling_strength=0.55, back_coupling=back)
        tt = n.cascade_from("thwaites")
        print(f"  {back:>16.2f} {f(tt['thwaites']):>12} {f(tt['ross']):>12}")

    print("\n  Thwaites tips at t=4 while Ross needs tau_slow=400 to respond,")
    print("  so no amount of Ross->Thwaites stabilisation arrives in time.")
    print("  A protective coupling is worthless if it is slower than the")
    print("  collapse it is meant to prevent — which is why the literature")
    print("  makes these interactions contingent on rate and delay, not just")
    print("  on sign.")

    print("\n" + "-" * 70)
    print("5. WHY IT MATTERS: consequence, not proximity")
    print("-" * 70)
    net = thwaites_ross_network(coupling_strength=0.55)
    t, x, h = net.simulate(
        lambda tt: np.array([0.45 if tt < 100 else 0.0, 0.0]),
        years=4000, dt=0.5)
    total = sum(e.consequence for e in net.elements)
    for when in (0, 200, 1000, 4000):
        i = min(int(when / 0.5), len(t) - 1)
        realised = net.consequence_realised(h[i])
        print(f"  t={t[i]:>6.0f}  x={np.array2string(x[i], precision=2):<16} "
              f"consequence realised {realised:>6.2f} / {total:.2f}")
    print("\n  Thwaites commits early and cheaply. Ross carries 95% of the "
          "cost\n  and arrives last — long after the trigger is gone.")
