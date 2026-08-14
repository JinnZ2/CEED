# Tests/test_universal_model.py

import numpy as np

from CEED_universal_model import (
    CEEDSystem,
    FeedbackLoop,
    SystemState,
    create_climate_system,
    create_financial_system,
)


def test_positive_feedback_amplifies():
    loop = FeedbackLoop(name="test", polarity="positive", strength=0.5)
    assert loop.apply(10.0) > 10.0


def test_negative_feedback_damps():
    loop = FeedbackLoop(name="test", polarity="negative", strength=0.5)
    assert loop.apply(10.0) < 10.0


def test_feedback_saturates_above_threshold():
    """Past the saturation threshold, effective strength is reduced."""
    loop = FeedbackLoop(
        name="test", polarity="positive", strength=0.5, saturation_threshold=5.0
    )
    below = loop.apply(1.0) / 1.0  # growth multiplier below threshold
    above = loop.apply(9.0) / 9.0  # growth multiplier above threshold
    assert above < below


def test_stability_metric():
    state = SystemState(energy=0.0, retention=2.0, dissipation=1.0, buffer_capacity=1.0)
    assert state.stability_metric() == 2.0

    stable = SystemState(energy=0.0, retention=0.5, dissipation=1.0, buffer_capacity=1.0)
    assert stable.stability_metric() < 1.0


def test_stability_metric_zero_dissipation_is_infinite():
    state = SystemState(energy=1.0, retention=1.0, dissipation=0.0, buffer_capacity=1.0)
    assert state.stability_metric() == float("inf")


def test_state_classification_boundaries():
    system = CEEDSystem("test")

    system.state.energy = 0.0
    assert system.classify_state() == "stable"

    system.state.energy = system.warning_threshold + 1
    assert system.classify_state() == "stressed"

    system.state.energy = system.critical_threshold + 1
    assert system.classify_state() == "critical"

    system.state.energy = system.tipping_point + 1
    assert system.classify_state() == "tipping"


def test_retention_collapses_with_energy():
    """Core CEED premise: high-energy states are harder to hold."""
    system = CEEDSystem("test")
    assert system.compute_retention(50.0) < system.compute_retention(0.0)


def test_dissipation_grows_with_energy():
    system = CEEDSystem("test")
    assert system.compute_dissipation(50.0) > system.compute_dissipation(10.0)


def test_simulate_returns_full_history():
    system = create_climate_system()
    history = system.simulate(steps=25, external_forcing_fn=lambda t: 2.0)

    assert len(history) == 25
    assert [h["timestep"] for h in history] == list(range(25))
    assert all(np.isfinite(h["energy"]) for h in history)


def test_diagnose_reports_expected_keys():
    system = create_financial_system()
    report = system.diagnose()

    for key in ("name", "state", "energy", "stability", "buffer_remaining",
                "distance_to_tipping", "retention", "dissipation", "num_feedbacks"):
        assert key in report

    assert report["name"] == "Finance"
    assert report["num_feedbacks"] == 4


def test_buffer_depletes_only_past_warning_threshold():
    system = CEEDSystem("test")

    # Below the warning threshold the buffer is untouched.
    system.state.energy = 1.0
    system.update(external_forcing=0.0)
    assert system.state.buffer_capacity == 1.0

    # Above it, the buffer erodes.
    system.state.energy = system.warning_threshold * 2
    system.update(external_forcing=0.0)
    assert system.state.buffer_capacity < 1.0
