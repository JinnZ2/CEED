"""Tests for the CEED convergence model."""

import pytest
import numpy as np

from simulation.convergence_model import (
    AnthropogenicForcing,
    ConvergencePredictor,
    ExtendedConvergencePredictor,
    SystemParameters,
    SYSTEMS,
)


# ── ConvergencePredictor (baseline) ───────────────────────────────

def test_model_initialization():
    predictor = ConvergencePredictor()
    assert isinstance(predictor.params.E_initial, dict)
    assert all(k in predictor.params.E_initial for k in SYSTEMS)


def test_energy_prediction_shape():
    predictor = ConvergencePredictor()
    t, solution = predictor.predict_convergence(years=1)
    assert len(t) == solution.shape[0]
    assert solution.shape[1] == 4


def test_phase_classification():
    predictor = ConvergencePredictor()
    t, solution = predictor.predict_convergence(years=1)
    total, phases = predictor.classify_phases(solution)
    assert len(phases) == len(total)
    # Phase 0 is the at-or-below-baseline band. It exists because thresholds
    # are now multiples of the baseline total rather than absolute values.
    assert all(0 <= p <= 4 for p in phases)


def test_energy_stays_finite():
    """ODE should not blow up over a long horizon."""
    predictor = ConvergencePredictor()
    t, solution = predictor.predict_convergence(years=10)
    assert np.all(np.isfinite(solution))


def test_dissipation_bounds_growth():
    """With no source and no retention, energy should decay toward zero."""
    params = SystemParameters()
    params.alpha_retention = {s: 0.0 for s in SYSTEMS}
    params.coupling = {}
    predictor = ConvergencePredictor(params)
    predictor.source_rate = lambda system, t: 0.0

    t, solution = predictor.predict_convergence(years=5)
    total_energy = np.sum(solution, axis=1)
    assert total_energy[-1] < total_energy[0]


def test_derivative_is_deterministic():
    """ODE RHS must return the same value for the same inputs."""
    predictor = ConvergencePredictor()
    E = [180.0, 92.5, 118.0, 110.0]
    d1 = predictor.energy_derivative(E, 0.5)
    d2 = predictor.energy_derivative(E, 0.5)
    assert d1 == d2


def test_coupling_is_conservative():
    """Energy transferred out of j must >= energy received by i (eta <= 1)."""
    params = SystemParameters()
    for (sys_i, sys_j), (c_ij, eta_ij) in params.coupling.items():
        assert 0 < eta_ij <= 1.0, (
            f"eta {sys_j}->{sys_i} = {eta_ij}, must be in (0, 1]")
        assert c_ij > 0, f"c {sys_j}->{sys_i} must be positive"


def test_coupling_produces_waste_heat():
    """With only coupling active, total energy must decrease (eta < 1)."""
    params = SystemParameters()
    params.alpha_retention = {s: 0.0 for s in SYSTEMS}
    params.lambda_dissipation = {s: 0.0 for s in SYSTEMS}
    params.gamma_nonlinear = {s: 0.0 for s in SYSTEMS}
    predictor = ConvergencePredictor(params)
    predictor.source_rate = lambda system, t: 0.0

    E = [180.0, 92.5, 118.0, 110.0]
    dE = predictor.energy_derivative(E, 0.0)
    assert sum(dE) < 0, f"Total dE/dt = {sum(dE):.6f}, should be < 0"


# ── ExtendedConvergencePredictor ──────────────────────────────────

def test_extended_returns_events():
    predictor = ExtendedConvergencePredictor()
    t, solution, events = predictor.predict_convergence(years=1)
    assert isinstance(events, list)
    assert len(t) == solution.shape[0]


def test_extended_no_external_matches_base():
    """Without external events, extended model should match baseline."""
    base = ConvergencePredictor()
    ext = ExtendedConvergencePredictor()

    t_b, sol_b = base.predict_convergence(years=1)
    t_e, sol_e, events = ext.predict_convergence(years=1, include_external=False)

    np.testing.assert_allclose(sol_b, sol_e, rtol=1e-10)


def test_extended_stays_finite():
    predictor = ExtendedConvergencePredictor()
    np.random.seed(42)
    t, solution, _ = predictor.predict_convergence(years=5)
    assert np.all(np.isfinite(solution))


# ── AnthropogenicForcing ─────────────────────────────────────────

def test_anthropogenic_forcing_grows():
    """Release rate should increase over time."""
    af = AnthropogenicForcing()
    assert af.release_rate(10) > af.release_rate(0)
    assert af.release_rate(20) > af.release_rate(10)


def test_anthropogenic_forcing_peaks():
    """With peak_year set, growth should saturate (logistic)."""
    af = AnthropogenicForcing(peak_year=50.0)
    # Growth rate should decelerate: the jump from t=0->50 should be
    # larger than the jump from t=200->250 (saturation)
    delta_early = af.release_rate(50) - af.release_rate(0)
    delta_late = af.release_rate(250) - af.release_rate(200)
    assert delta_late < delta_early * 0.1  # late growth is <10% of early


def test_anthropogenic_forcing_no_peak():
    """Without peak_year, growth is pure exponential."""
    af = AnthropogenicForcing(peak_year=None)
    expected = af.A_0 * np.exp(af.growth_rate * 10)
    assert abs(af.release_rate(10) - expected) < 1e-10


def test_anthropogenic_cumulative_positive():
    """Cumulative release must be positive and increasing."""
    af = AnthropogenicForcing()
    assert af.cumulative_release(1) > 0
    assert af.cumulative_release(10) > af.cumulative_release(5)


# ── State-dependent parameters ────────────────────────────────────

def test_no_anthro_returns_base_params():
    """Without anthropogenic forcing, effective params = baseline."""
    params = SystemParameters()  # anthro is None
    predictor = ConvergencePredictor(params)
    E = [180.0, 92.5, 118.0, 110.0]

    for s in SYSTEMS:
        assert predictor.effective_alpha(s, E, 0.0) == params.alpha_retention[s]

    for key, (c_base, eta_base) in params.coupling.items():
        c_eff, eta_eff = predictor.effective_coupling(key[0], key[1], E, 0.0)
        assert c_eff == c_base
        assert eta_eff == eta_base


def test_anthro_increases_atm_retention():
    """Anthropogenic forcing should increase atmospheric retention."""
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    predictor = ConvergencePredictor(params)
    E = [180.0, 92.5, 118.0, 110.0]

    alpha_base = params.alpha_retention['atmospheric']
    alpha_eff = predictor.effective_alpha('atmospheric', E, 0.0)
    assert alpha_eff > alpha_base


def test_anthro_increases_coupling_with_gradient():
    """Larger energy gradient should strengthen coupling."""
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    predictor = ConvergencePredictor(params)

    E_small_gradient = [100.0, 100.0, 100.0, 100.0]
    E_large_gradient = [200.0, 50.0, 150.0, 100.0]

    c_small, _ = predictor.effective_coupling(
        'atmospheric', 'solar', E_small_gradient, 0.0)
    c_large, _ = predictor.effective_coupling(
        'atmospheric', 'solar', E_large_gradient, 0.0)

    assert c_large > c_small


def test_anthro_efficiency_clamped():
    """Conversion efficiency must stay in [0.01, 0.95]."""
    params = SystemParameters(
        anthropogenic=AnthropogenicForcing(efficiency_shift=1.0))
    predictor = ConvergencePredictor(params)

    # Very high energy should push efficiency toward cap
    E_high = [1000.0, 1000.0, 1000.0, 1000.0]
    for (sys_i, sys_j) in params.coupling:
        _, eta = predictor.effective_coupling(sys_i, sys_j, E_high, 0.0)
        assert 0.01 <= eta <= 0.95


def test_anthro_model_stays_finite():
    """ODE with anthropogenic forcing should not blow up."""
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    predictor = ConvergencePredictor(params)
    t, solution = predictor.predict_convergence(years=20)
    assert np.all(np.isfinite(solution))


def test_anthro_increases_total_energy():
    """Anthropogenic forcing should result in higher total energy than baseline."""
    base = ConvergencePredictor(SystemParameters())
    anthro = ConvergencePredictor(
        SystemParameters(anthropogenic=AnthropogenicForcing()))

    _, sol_b = base.predict_convergence(years=10)
    _, sol_a = anthro.predict_convergence(years=10)

    E_base_final = np.sum(sol_b[-1])
    E_anthro_final = np.sum(sol_a[-1])
    assert E_anthro_final > E_base_final


def test_long_horizon_stability():
    """Model must stay finite and non-negative over 200 years."""
    for params in [SystemParameters(),
                   SystemParameters(anthropogenic=AnthropogenicForcing())]:
        predictor = ConvergencePredictor(params)
        t, solution = predictor.predict_convergence(years=200)
        assert np.all(np.isfinite(solution)), "Solution went infinite"
        # Energy should be non-negative (soft floor keeps it near zero)
        assert np.all(solution > -1.0), (
            f"Energy went significantly negative: min={solution.min():.1f}"
        )
