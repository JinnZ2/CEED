# Tests/test_universal_model.py

import numpy as np
import pytest

from CEED_universal_model import (
    CEEDSystem,
    FeedbackLoop,
    SystemState,
    create_climate_system,
    create_financial_system,
)


# --- feedback loops --------------------------------------------------------

def test_positive_feedback_adds_energy():
    loop = FeedbackLoop(name="test", polarity="positive", strength=0.5)
    assert loop.rate(10.0) > 0


def test_negative_feedback_removes_energy():
    loop = FeedbackLoop(name="test", polarity="negative", strength=0.5)
    assert loop.rate(10.0) < 0


def test_unsaturated_feedback_is_linear_in_energy():
    loop = FeedbackLoop(name="test", polarity="positive", strength=0.4)
    assert loop.rate(20.0) == pytest.approx(2.0 * loop.rate(10.0))
    assert loop.rate(10.0) == pytest.approx(0.4 * 10.0)


def test_feedback_vanishes_at_zero_energy():
    loop = FeedbackLoop(name="test", polarity="positive", strength=0.5)
    assert loop.rate(0.0) == pytest.approx(0.0)


def test_saturation_suppresses_feedback_at_high_energy():
    """sigma = 1/(1+(E/E_sat)^2), so the per-unit gain falls past E_sat."""
    loop = FeedbackLoop(name="test", polarity="positive", strength=0.5,
                        saturation_energy=5.0)
    gain_low = loop.rate(1.0) / 1.0
    gain_high = loop.rate(50.0) / 50.0
    assert gain_high < gain_low


def test_saturation_at_the_scale_energy_halves_the_gain():
    """At E = E_sat, sigma = 1/(1+1) = 0.5 exactly."""
    loop = FeedbackLoop(name="test", polarity="positive", strength=1.0,
                        saturation_energy=8.0)
    assert loop.rate(8.0) == pytest.approx(0.5 * 1.0 * 8.0)


def test_saturating_feedback_rate_eventually_decays():
    """s*E/(1+(E/Es)^2) -> 0 as E grows without bound."""
    loop = FeedbackLoop(name="test", polarity="positive", strength=1.0,
                        saturation_energy=5.0)
    assert loop.rate(1e6) < loop.rate(5.0)


# --- dissipation -----------------------------------------------------------

def test_dissipation_grows_with_energy():
    system = CEEDSystem("test")
    assert system.dissipation_rate(50.0) > system.dissipation_rate(10.0)


def test_dissipation_vanishes_at_zero():
    assert CEEDSystem("test").dissipation_rate(0.0) == pytest.approx(0.0)


def test_dissipation_matches_its_formula():
    system = CEEDSystem("test")
    E = 40.0
    expected = (system.dissipation_linear * E
                + system.dissipation_nonlinear * E ** system.dissipation_exponent)
    assert system.dissipation_rate(E) == pytest.approx(expected)


def test_dissipation_is_superlinear():
    """Doubling energy more than doubles dissipation, thanks to the E^1.5 term."""
    system = CEEDSystem("test")
    assert system.dissipation_rate(80.0) > 2.0 * system.dissipation_rate(40.0)


# --- state and derivative --------------------------------------------------

def test_dE_dt_combines_forcing_feedback_and_dissipation():
    system = CEEDSystem("test")
    system.add_feedback(FeedbackLoop("f", "positive", 0.2))
    E = 30.0
    expected = 5.0 + system.total_feedback_rate(E) - system.dissipation_rate(E)
    assert system.dE_dt(E, external_forcing=5.0) == pytest.approx(expected)


def test_stability_metric_above_one_when_feedback_is_net_positive():
    system = CEEDSystem("test")
    system.add_feedback(FeedbackLoop("amp", "positive", 0.5))
    system.state.energy = 10.0
    assert system.stability_metric() > 1.0


def test_stability_metric_below_one_when_feedback_is_net_negative():
    system = CEEDSystem("test")
    system.add_feedback(FeedbackLoop("damp", "negative", 0.5))
    system.state.energy = 10.0
    assert system.stability_metric() < 1.0


def test_stability_metric_is_one_with_no_feedbacks():
    system = CEEDSystem("test")
    system.state.energy = 10.0
    assert system.stability_metric() == pytest.approx(1.0)


def test_stability_metric_infinite_when_dissipation_zero():
    system = CEEDSystem("test")
    system.state.energy = 0.0  # D(0) = 0
    assert system.stability_metric() == float('inf')


def test_system_state_defaults_to_full_buffer():
    assert SystemState(energy=1.0).buffer_capacity == 1.0


# --- classification --------------------------------------------------------

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


def test_all_four_classifications_are_reachable():
    """Guards the defect where thresholds sat outside the reachable range."""
    system = CEEDSystem("test")
    seen = set()
    for E in (0.0, 150.0, 250.0, 400.0):
        system.state.energy = E
        seen.add(system.classify_state())
    assert seen == {"stable", "stressed", "critical", "tipping"}


# --- buffer ----------------------------------------------------------------

def test_buffer_intact_below_warning_threshold():
    system = CEEDSystem("test")
    system.state.energy = 1.0
    system.update(external_forcing=0.0)
    assert system.state.buffer_capacity == 1.0


def test_buffer_depletes_above_warning_threshold():
    system = CEEDSystem("test")
    system.state.energy = system.warning_threshold * 2
    system.update(external_forcing=0.0)
    assert system.state.buffer_capacity < 1.0


def test_buffer_never_goes_negative():
    system = CEEDSystem("test")
    system.state.energy = 10_000.0
    for _ in range(500):
        system.update(external_forcing=0.0, dt=0.5)
    assert system.state.buffer_capacity >= 0.0


# --- simulation ------------------------------------------------------------

def test_simulate_returns_full_history():
    system = create_climate_system()
    history = system.simulate(steps=25, external_forcing_fn=lambda t: 2.0)

    assert len(history) == 25
    assert [h["timestep"] for h in history] == list(range(25))
    assert all(np.isfinite(h["energy"]) for h in history)


def test_simulate_records_time_consistent_with_dt():
    system = CEEDSystem("test")
    history = system.simulate(steps=5, dt=0.25)
    assert [h["time"] for h in history] == pytest.approx([0.0, 0.25, 0.5, 0.75, 1.0])


def test_zero_forcing_from_zero_energy_stays_put():
    """E=0 is a fixed point: no forcing, no feedback, no dissipation."""
    system = create_climate_system()
    history = system.simulate(steps=10)
    assert history[-1]["energy"] == pytest.approx(0.0)


def test_diagnose_reports_expected_keys():
    system = create_financial_system()
    report = system.diagnose()

    for key in ("name", "state", "energy", "stability", "buffer_remaining",
                "distance_to_tipping", "feedback_rate", "dissipation_rate",
                "num_feedbacks"):
        assert key in report

    assert report["name"] == "Finance"
    assert report["num_feedbacks"] == 4


def test_example_systems_carry_both_polarities():
    for build in (create_climate_system, create_financial_system):
        polarities = {fb.polarity for fb in build().feedbacks}
        assert polarities == {"positive", "negative"}
