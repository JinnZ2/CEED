"""Tests for the CEED convergence model."""

import pytest
import numpy as np

from simulation.convergence_model import (
    ConvergencePredictor,
    ExtendedConvergencePredictor,
    SystemParameters,
    SYSTEMS,
)


# ── ConvergencePredictor ──────────────────────────────────────────

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
    assert all(1 <= p <= 4 for p in phases)


def test_energy_stays_finite():
    """ODE should not blow up over a long horizon."""
    predictor = ConvergencePredictor()
    t, solution = predictor.predict_convergence(years=10)
    assert np.all(np.isfinite(solution))


def test_dissipation_bounds_growth():
    """With no source, energy should decay toward zero."""
    params = SystemParameters()
    # Zero out all sources by setting initial energies only
    predictor = ConvergencePredictor(params)

    # Override source_rate to return 0
    predictor.source_rate = lambda system, t: 0.0

    # Set coupling to zero as well
    params.coupling = {}
    predictor.params = params

    # Start with some energy
    t, solution = predictor.predict_convergence(years=5)
    total_energy = np.sum(solution, axis=1)
    # Energy should decrease (dissipation only, no source)
    assert total_energy[-1] < total_energy[0]


def test_derivative_is_deterministic():
    """ODE RHS must return the same value for the same inputs."""
    predictor = ConvergencePredictor()
    E = [180.0, 92.5, 118.0, 110.0]
    d1 = predictor.energy_derivative(E, 0.5)
    d2 = predictor.energy_derivative(E, 0.5)
    assert d1 == d2


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
