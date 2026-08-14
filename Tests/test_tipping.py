# Tests/test_tipping.py

"""
Tests for the tipping-element primitive.

The three that matter are `test_stays_tipped_after_forcing_is_removed`,
`test_hysteresis_loop_has_width`, and
`test_commitment_precedes_visible_change` — they assert the properties the
convergence model provably lacks. `test_convergence_model_has_one_attractor`
pins that contrast so it cannot drift silently.
"""

import numpy as np
import pytest

from simulation.tipping import (
    FOLD_F,
    FOLD_X,
    TippingElement,
    hysteresis_loop,
    recovery_rate,
    rolling_autocorrelation,
    rolling_variance,
)


# --- structure -------------------------------------------------------------

def test_fold_constants_match_the_analytic_values():
    """dx/dt = x - x^3 + F folds at x = +/-1/sqrt(3), F = -/+2/(3 sqrt 3)."""
    assert FOLD_X == pytest.approx(1.0 / np.sqrt(3.0))
    assert FOLD_F == pytest.approx(2.0 / (3.0 * np.sqrt(3.0)))
    assert FOLD_F == pytest.approx(0.3849, abs=1e-4)


def test_two_stable_branches_below_the_fold():
    stable = TippingElement.stable_equilibria(0.0)
    assert len(stable) == 2
    assert stable[0] == pytest.approx(-1.0)
    assert stable[1] == pytest.approx(1.0)


def test_one_branch_beyond_the_fold():
    assert len(TippingElement.stable_equilibria(0.6)) == 1
    assert len(TippingElement.stable_equilibria(-0.6)) == 1


def test_bistability_window_matches_the_fold():
    assert TippingElement.is_bistable(0.0)
    assert TippingElement.is_bistable(FOLD_F * 0.99)
    assert not TippingElement.is_bistable(FOLD_F * 1.01)
    assert not TippingElement.is_bistable(-FOLD_F * 1.01)


def test_equilibria_actually_solve_the_equation():
    for F in (-0.3, 0.0, 0.2, 0.5):
        for x in TippingElement.equilibria(F):
            assert (x - x ** 3 + F) == pytest.approx(0.0, abs=1e-9)


# --- the properties the convergence model lacks ----------------------------

def test_stays_tipped_after_forcing_is_removed():
    """
    Hysteresis. Push past the fold, remove the forcing entirely, integrate
    far longer than either timescale. The system must NOT return.
    """
    el = TippingElement(tau_fast=1.0, tau_slow=50.0, x0=-1.0)
    t, x, h = el.simulate(lambda tt: 0.45 if tt < 100 else 0.0,
                          years=2000, dt=0.1)

    assert x[0] == pytest.approx(-1.0)
    assert x[-1] == pytest.approx(1.0, abs=1e-3)
    assert h[-1] == pytest.approx(1.0, abs=1e-3)


def test_sub_threshold_forcing_does_not_tip():
    """The mirror: below the fold the system returns to where it started."""
    el = TippingElement(tau_fast=1.0, tau_slow=50.0, x0=-1.0)
    t, x, h = el.simulate(lambda tt: 0.30 if tt < 100 else 0.0,
                          years=2000, dt=0.1)
    assert x[-1] < 0
    assert x[-1] == pytest.approx(-1.0, abs=1e-3)


def test_hysteresis_loop_has_width():
    """Up-ramp and down-ramp must switch at different forcings."""
    F_up, x_up, F_down, x_down = hysteresis_loop()

    up_switch = F_up[np.argmax(x_up > 0)]
    down_switch = F_down[np.argmax(x_down < 0)]

    assert up_switch > 0
    assert down_switch < 0
    assert up_switch - down_switch > 0.5


def test_switch_points_are_near_the_analytic_folds():
    """A slow ramp should switch close to +/-FOLD_F, slightly late."""
    F_up, x_up, F_down, x_down = hysteresis_loop(ramp_years=8000.0, dt=0.5)
    up_switch = F_up[np.argmax(x_up > 0)]
    assert up_switch == pytest.approx(FOLD_F, abs=0.03)
    assert up_switch >= FOLD_F  # finite ramp rate lags the static fold


def test_commitment_precedes_visible_change():
    """
    The load-bearing property. The internal state must cross while the
    observable has barely moved — the system is decided before it looks it.
    """
    el = TippingElement(tau_fast=1.0, tau_slow=50.0, x0=-1.0)
    t, x, h = el.simulate(lambda tt: 0.45 if tt < 100 else 0.0,
                          years=1200, dt=0.1)

    i_commit = np.where(x > 0)[0][0]
    progress = (h[i_commit] - h[0]) / (h.max() - h[0])

    assert progress < 0.10          # under 10% of the eventual change
    lag = el.commitment_lag(t, x, h)
    assert lag is not None and lag > 10.0


def test_commitment_lag_scales_with_timescale_separation():
    """Slower response -> longer blind window."""
    lags = []
    for tau_slow in (20.0, 100.0):
        el = TippingElement(tau_fast=1.0, tau_slow=tau_slow, x0=-1.0)
        t, x, h = el.simulate(lambda tt: 0.45 if tt < 100 else 0.0,
                              years=3000, dt=0.1)
        lags.append(el.commitment_lag(t, x, h))
    assert lags[1] > lags[0]


def test_commitment_lag_is_none_when_nothing_tips():
    el = TippingElement(tau_fast=1.0, tau_slow=50.0, x0=-1.0)
    t, x, h = el.simulate(lambda tt: 0.0, years=500, dt=0.1)
    assert el.commitment_lag(t, x, h) is None


def test_is_committed_tracks_the_unstable_root():
    el = TippingElement(x0=-1.0)
    assert not el.is_committed(0.0, x=-0.5)
    assert el.is_committed(0.0, x=+0.5)
    assert el.is_committed(0.5)          # past the fold, only one branch
    assert not el.is_committed(-0.5)


# --- early warning signals -------------------------------------------------

def test_recovery_slows_approaching_the_fold():
    """Critical slowing down: the recovery rate goes to zero at the fold."""
    el = TippingElement()
    rates = [recovery_rate(el, F) for F in (0.0, 0.2, 0.3, 0.37, 0.384)]
    assert all(a > b for a, b in zip(rates, rates[1:]))
    assert rates[-1] < 0.15


def test_recovery_rate_rejects_a_missing_branch():
    el = TippingElement()
    with pytest.raises(ValueError):
        recovery_rate(el, 0.0, branch='nonsense')


def test_autocorrelation_rises_before_tipping():
    """
    The standard early-warning signal, on a noisy trajectory. Noise is
    pre-generated, never drawn inside the derivative.
    """
    rng = np.random.default_rng(20260814)
    el = TippingElement(tau_fast=1.0, tau_slow=50.0, x0=-1.0)
    steps = int(3000 / 0.1) + 1
    noise = rng.normal(0.0, 0.05, steps)

    # Slow ramp from well below the fold up to just under it.
    t, x, h = el.simulate(lambda tt: 0.38 * (tt / 3000.0),
                          years=3000, dt=0.1, noise=noise)

    ac = rolling_autocorrelation(x, window=400)
    early = np.nanmean(ac[2000:6000])
    late = np.nanmean(ac[-6000:-2000])
    assert late > early


def test_rolling_statistics_shapes_and_guards():
    series = np.arange(50, dtype=float)
    assert len(rolling_autocorrelation(series, 10)) == len(series)
    assert len(rolling_variance(series, 10)) == len(series)
    assert np.isnan(rolling_autocorrelation(series, 10)[:10]).all()
    for fn in (rolling_autocorrelation, rolling_variance):
        with pytest.raises(ValueError):
            fn(series, 2)


def test_rolling_variance_detects_a_variance_increase():
    rng = np.random.default_rng(7)
    calm = rng.normal(0, 0.1, 500)
    wild = rng.normal(0, 1.0, 500)
    v = rolling_variance(np.concatenate([calm, wild]), window=100)
    assert np.nanmean(v[900:]) > 10 * np.nanmean(v[300:450])


# --- construction and determinism ------------------------------------------

def test_rejects_inverted_timescales():
    """Without tau_slow > tau_fast there is no commitment window."""
    with pytest.raises(ValueError):
        TippingElement(tau_fast=10.0, tau_slow=1.0)
    with pytest.raises(ValueError):
        TippingElement(tau_fast=-1.0, tau_slow=50.0)


def test_derivative_is_deterministic():
    """Project rule: no random() inside any RHS."""
    el = TippingElement()
    a = el.derivative([0.3, -0.2], F=0.1, noise=0.02)
    b = el.derivative([0.3, -0.2], F=0.1, noise=0.02)
    assert np.array_equal(a, b)


def test_simulate_rejects_short_noise_array():
    el = TippingElement()
    with pytest.raises(ValueError):
        el.simulate(lambda t: 0.0, years=100, dt=0.1, noise=np.zeros(5))


def test_trajectory_is_finite_and_correctly_shaped():
    el = TippingElement()
    t, x, h = el.simulate(lambda tt: 0.2, years=100, dt=0.1)
    assert len(t) == len(x) == len(h) == 1001
    assert np.all(np.isfinite(x)) and np.all(np.isfinite(h))


# --- the contrast being asserted -------------------------------------------

def test_convergence_model_has_one_attractor():
    """
    Pins the property that motivates this module. If the convergence model
    ever gains genuine bistability, this test fails and the docs claiming
    it has none need updating with it.
    """
    from simulation.convergence_model import (SYSTEMS, ConvergencePredictor,
                                              SystemParameters)

    base = SystemParameters().E_initial
    finals = []
    for mult in (0.5, 1.0, 1.6, 3.0):
        pr = ConvergencePredictor(SystemParameters(
            E_initial={s: mult * base[s] for s in SYSTEMS}))
        _, sol = pr.predict_convergence(years=200)
        finals.append(float(np.sum(sol[-1])))

    assert max(finals) - min(finals) < 1.0
