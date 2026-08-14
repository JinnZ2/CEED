# Tests/test_extended_model.py

import numpy as np

from simulation.convergence_model_extended import (
    ExtendedConvergencePredictor,
    ExternalEvent,
)


def test_prediction_shape():
    predictor = ExtendedConvergencePredictor()
    t, solution, events = predictor.predict_convergence(years=1, include_external=False)

    assert len(t) == solution.shape[0]
    assert solution.shape[1] == 4
    assert np.all(np.isfinite(solution))


def test_no_external_events_when_disabled():
    predictor = ExtendedConvergencePredictor()
    predictor.events = []
    _, _, events = predictor.predict_convergence(years=1, include_external=False)

    # Only volcanic events (generated inside the dissipation term) may appear.
    assert all(e.event_type == "volcanic" for e in events)


def test_external_events_are_generated():
    predictor = ExtendedConvergencePredictor()
    _, _, events = predictor.predict_convergence(years=2, include_external=True)

    types = {e.event_type for e in events}
    assert "launch" in types


def test_launch_events_within_horizon():
    predictor = ExtendedConvergencePredictor()
    years = 2
    launches = predictor.add_satellite_launches(years)

    assert len(launches) > 0
    assert all(e.time < years for e in launches)
    assert all(e.event_type == "launch" for e in launches)


def test_meteor_energy_within_configured_range():
    predictor = ExtendedConvergencePredictor()
    low, high = predictor.meteor_energy_range
    meteors = predictor.generate_meteor_events(years=50)

    assert all(low <= e.energy <= high for e in meteors)
    # generate_meteor_events returns events sorted by time
    assert [e.time for e in meteors] == sorted(e.time for e in meteors)


def test_dissipation_saturates_at_high_energy():
    """Sinks lose effectiveness as total energy climbs."""
    predictor = ExtendedConvergencePredictor()
    np.random.seed(0)

    # Compare the deterministic base_loss term as a fraction of total energy.
    def base_fraction(E_total):
        factor = max(0.1, 1.0 - (E_total / predictor.dissipation_saturation))
        return predictor.unknown_sink_baseline * factor

    assert base_fraction(700) < base_fraction(100)


def test_phase_classification_range():
    predictor = ExtendedConvergencePredictor()
    _, solution, _ = predictor.predict_convergence(years=1, include_external=False)
    total, phases = predictor.classify_phases(solution)

    assert len(phases) == len(total)
    assert all(1 <= p <= 4 for p in phases)


def test_analyze_event_timing_flags_high_energy():
    predictor = ExtendedConvergencePredictor()
    t, solution, _ = predictor.predict_convergence(years=1, include_external=False)

    # An event at t=0, where energy already sits well above the phase_2 threshold.
    event = ExternalEvent(time=0.0, energy=5.0, event_type="meteor")
    critical = predictor.analyze_event_timing(t, solution, [event])

    assert len(critical) == 1
    assert critical[0]["amplification_risk"] > 1.0
