# Tests/test_extended_model.py

import numpy as np
import pytest

from simulation.convergence_model import (
    SYSTEMS,
    AnthropogenicForcing,
    ConvergencePredictor,
    ExtendedConvergencePredictor,
    ExternalEvent,
    SystemParameters,
)


# --- integration shape -----------------------------------------------------

def test_prediction_shape():
    predictor = ExtendedConvergencePredictor()
    t, solution, events = predictor.predict_convergence(years=1, include_external=False)

    assert len(t) == solution.shape[0]
    assert solution.shape[1] == len(SYSTEMS)
    assert np.all(np.isfinite(solution))


def test_no_events_when_disabled():
    predictor = ExtendedConvergencePredictor()
    _, _, events = predictor.predict_convergence(years=1, include_external=False)
    assert events == []


def test_external_events_are_generated():
    predictor = ExtendedConvergencePredictor()
    _, _, events = predictor.predict_convergence(years=2, include_external=True)
    assert "launch" in {e.event_type for e in events}


def test_events_are_pregenerated_so_the_rhs_is_deterministic():
    """
    Events are generated before integration, so evaluating the derivative
    twice at the same state must give the same answer. A stochastic RHS
    would silently break the ODE solver's error control.
    """
    predictor = ExtendedConvergencePredictor()
    predictor.predict_convergence(years=1, include_external=True)

    E0 = [predictor.params.E_initial[s] for s in SYSTEMS]
    first = predictor.energy_derivative(E0, 0.5)
    second = predictor.energy_derivative(E0, 0.5)
    assert first == pytest.approx(second)


# --- event generation ------------------------------------------------------

def test_launch_events_within_horizon():
    predictor = ExtendedConvergencePredictor()
    years = 2
    events = predictor._generate_events(years)
    launches = [e for e in events if e.event_type == "launch"]

    assert len(launches) > 0
    assert all(0.0 <= e.time <= years for e in launches)


def test_meteor_energy_within_configured_range():
    predictor = ExtendedConvergencePredictor()
    low, high = predictor.meteor_energy_range
    meteors = [e for e in predictor._generate_events(50) if e.event_type == "meteor"]

    assert len(meteors) > 0
    assert all(low <= e.energy <= high for e in meteors)


def test_generated_events_are_time_ordered():
    predictor = ExtendedConvergencePredictor()
    events = predictor._generate_events(10)
    times = [e.time for e in events]
    assert times == sorted(times)


# --- phase classification --------------------------------------------------

def test_phase_thresholds_are_relative_to_baseline():
    """
    Thresholds must scale with the baseline total, not be absolute. The
    earlier absolute values (120/150/200/300) sat below a baseline of 500.5,
    so every run classified as phase 4 from t=0.
    """
    params = SystemParameters()
    base = params.E_baseline
    assert base == pytest.approx(500.5)

    thresholds = params.phase_thresholds
    assert thresholds['phase_1'] == pytest.approx(1.20 * base)
    assert thresholds['phase_4'] == pytest.approx(3.00 * base)
    assert thresholds['phase_1'] > base


def test_baseline_run_is_not_pinned_at_phase_four():
    """The regression guard for the unreachable-threshold defect."""
    predictor = ConvergencePredictor()
    _, solution = predictor.predict_convergence(years=3)
    _, phases = predictor.classify_phases(solution)
    assert max(phases) < 4


def test_phase_classification_range():
    predictor = ExtendedConvergencePredictor()
    _, solution, _ = predictor.predict_convergence(years=1, include_external=False)
    total, phases = predictor.classify_phases(solution)

    assert len(phases) == len(total)
    assert all(0 <= p <= 4 for p in phases)


def test_phases_track_total_energy_monotonically():
    predictor = ConvergencePredictor()
    thresholds = predictor.params.phase_thresholds
    fake = np.array([
        [100.0, 100.0, 100.0, 100.0],                  # 400, below phase_1
        [thresholds['phase_1'] / 4] * 4,               # exactly phase_1
        [thresholds['phase_3'] / 4] * 4,               # phase_3
        [thresholds['phase_4'] / 4 + 10] * 4,          # above phase_4
    ])
    _, phases = predictor.classify_phases(fake)
    assert phases == [0, 1, 3, 4]


# --- event timing ----------------------------------------------------------

def test_event_timing_is_measured_against_baseline():
    """
    Amplification risk is E_at_event / E_baseline. Keyed to an absolute phase
    threshold this test was degenerate — it flagged either every event or none.
    """
    predictor = ExtendedConvergencePredictor()
    t, solution, _ = predictor.predict_convergence(years=1, include_external=False)

    event = ExternalEvent(time=0.0, energy=5.0, event_type="meteor")
    # threshold_ratio=0 flags everything, so risk is always reported
    critical = predictor.analyze_event_timing(t, solution, [event], threshold_ratio=0.0)

    assert len(critical) == 1
    expected = np.sum(solution, axis=1)[0] / predictor.params.E_baseline
    assert critical[0]['amplification_risk'] == pytest.approx(expected, rel=1e-6)


def test_quiet_run_flags_no_critical_events():
    """A run that never leaves the baseline band should flag nothing."""
    predictor = ExtendedConvergencePredictor()
    t, solution, events = predictor.predict_convergence(years=3, include_external=True)
    critical = predictor.analyze_event_timing(t, solution, events)

    total = np.sum(solution, axis=1)
    if total.max() < predictor.params.phase_thresholds['phase_1']:
        assert critical == []


def test_negative_energy_events_are_skipped():
    predictor = ExtendedConvergencePredictor()
    t, solution, _ = predictor.predict_convergence(years=1, include_external=False)
    sink = ExternalEvent(time=0.0, energy=-10.0, event_type="volcanic")
    assert predictor.analyze_event_timing(t, solution, [sink],
                                          threshold_ratio=0.0) == []


# --- anthropogenic forcing -------------------------------------------------

def test_anthropogenic_forcing_adds_energy():
    base = ConvergencePredictor(SystemParameters())
    anthro = ConvergencePredictor(
        SystemParameters(anthropogenic=AnthropogenicForcing()))

    _, sol_base = base.predict_convergence(years=10)
    _, sol_anthro = anthro.predict_convergence(years=10)

    assert np.sum(sol_anthro[-1]) > np.sum(sol_base[-1])


def test_efficiency_shift_uses_fractional_excess():
    """
    eta_eff = eta_base * (1 + delta*(E_tot - E_ref)/E_ref). Without the
    /E_ref the response was inflated by a factor of E_ref (~500).
    """
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    predictor = ConvergencePredictor(params)

    E_ref = params.E_baseline
    at_baseline = [params.E_initial[s] for s in SYSTEMS]
    _, eta_base = predictor.effective_coupling('oceanic', 'atmospheric',
                                               at_baseline, 0.0)

    # Double the total energy -> fractional excess of exactly 1.0
    doubled = [2.0 * e for e in at_baseline]
    _, eta_doubled = predictor.effective_coupling('oceanic', 'atmospheric',
                                                  doubled, 0.0)

    delta = params.anthropogenic.efficiency_shift
    assert eta_doubled == pytest.approx(eta_base * (1.0 + delta), rel=1e-6)
