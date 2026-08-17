# Tests/test_forecast.py

"""
Tests for the pre-registered forecast.

The load-bearing property is that the forecast is genuinely out of sample:
`test_forecast_uses_no_observation_after_initialisation` poisons every
post-2024 observation and asserts the predictions do not move. Everything
else about pre-registration is decoration if that fails.
"""

import numpy as np
import pytest

from simulation.convergence_model import SYSTEMS, AnthropogenicForcing, SystemParameters
from simulation.forecast import (
    INTERVAL_SIGMA,
    LAST_OBSERVED_YEAR,
    Prediction,
    SkillPrior,
    falsification_criteria,
    run_forecast,
    skill_priors,
    to_record,
    trend_comparison,
)
from simulation.unit_bridge import ALL_OBSERVED


@pytest.fixture(scope="module")
def priors():
    return skill_priors()


@pytest.fixture(scope="module")
def predictions(priors):
    return run_forecast(priors=priors)


# --- the out-of-sample guarantee -------------------------------------------

def test_forecast_uses_no_observation_after_initialisation(monkeypatch, priors):
    """
    Corrupt every observation after the initialisation year. A genuinely
    out-of-sample forecast is unchanged.
    """
    before = run_forecast(priors=priors)

    for system, table in ALL_OBSERVED.items():
        poisoned = dict(table)
        for year in list(poisoned):
            if year > LAST_OBSERVED_YEAR:
                poisoned[year] = poisoned[year] * 10.0
        monkeypatch.setitem(ALL_OBSERVED, system, poisoned)

    after = run_forecast(priors=priors)
    assert [p.value for p in after] == pytest.approx([p.value for p in before])


def test_forecast_does_depend_on_the_initialisation_year(monkeypatch, priors):
    """The mirror: changing the 2024 anchor must move the forecast."""
    before = run_forecast(priors=priors)

    table = dict(ALL_OBSERVED['oceanic'])
    table[LAST_OBSERVED_YEAR] = table[LAST_OBSERVED_YEAR] * 1.5
    monkeypatch.setitem(ALL_OBSERVED, 'oceanic', table)

    after = run_forecast(priors=priors)
    b = [p.value for p in before if p.system == 'oceanic']
    a = [p.value for p in after if p.system == 'oceanic']
    assert a != pytest.approx(b)


def test_forecast_rejects_a_year_without_observations():
    with pytest.raises(ValueError):
        run_forecast(start_year=1850)


# --- shape and banding -----------------------------------------------------

def test_one_prediction_per_system_per_year(predictions):
    horizon = 6
    assert len(predictions) == horizon * len(SYSTEMS)
    years = sorted({p.year for p in predictions})
    assert years == list(range(LAST_OBSERVED_YEAR + 1,
                               LAST_OBSERVED_YEAR + horizon + 1))


def test_bands_are_assigned_relative_to_issue_year(priors):
    preds = run_forecast(issue_year=2026, priors=priors)
    band = {p.year: p.band for p in preds}
    assert band[2025] == "resolvable now"
    assert band[2026] == "resolving"
    assert band[2027] == "prospective"


def test_every_prediction_is_finite(predictions):
    assert all(np.isfinite([p.value, p.lower, p.upper]).all()
               for p in predictions)


# --- intervals -------------------------------------------------------------

def test_intervals_come_from_held_out_test_rmse(predictions, priors):
    for p in predictions:
        half = INTERVAL_SIGMA * priors[p.system].test_rmse
        assert p.upper - p.value == pytest.approx(half)
        assert p.value - p.lower == pytest.approx(half)


def test_intervals_bracket_the_prediction(predictions):
    assert all(p.lower < p.value < p.upper for p in predictions)


def test_oceanic_interval_is_widest_relative_to_its_skill(priors):
    """The subsystem that failed hardest out of sample gets the widest band."""
    assert priors['oceanic'].test_rmse > priors['atmospheric'].test_rmse


# --- declared expectations -------------------------------------------------

def test_expectation_thresholds():
    assert SkillPrior('x', 0.6, 1.0, 'u').expectation == "expect skill"
    assert SkillPrior('x', 0.1, 1.0, 'u').expectation == "marginal"
    assert "NO skill" in SkillPrior('x', -1.0, 1.0, 'u').expectation


def test_solar_is_declared_skilful_and_oceanic_is_not(priors):
    """
    The meta-prediction: validation says solar generalises and oceanic does
    not. If this ordering ever flips, the forecast's premise changed.
    """
    assert priors['solar'].test_r2 > 0.3
    assert priors['oceanic'].test_r2 < 0.0


def test_predictions_carry_their_declared_skill(predictions, priors):
    for p in predictions:
        assert p.test_r2 == pytest.approx(priors[p.system].test_r2)
        assert p.expectation == priors[p.system].expectation


# --- trend comparison ------------------------------------------------------

def test_cyclic_systems_are_flagged_not_comparable(predictions):
    trends = trend_comparison(predictions)
    assert trends['solar']['comparable'] is False
    assert trends['magnetic']['comparable'] is False
    assert trends['atmospheric']['comparable'] is True
    assert trends['oceanic']['comparable'] is True


def test_oceanic_forecast_undershoots_the_observed_trend(predictions):
    """
    The specific predicted failure: the model continues the ocean record at a
    fraction of its observed rate. Recorded so the miss is a prediction, not
    a post-hoc explanation.
    """
    oce = trend_comparison(predictions)['oceanic']
    assert 0.0 < oce['ratio'] < 0.5
    assert oce['forecast_per_year'] < oce['observed_per_year']


def test_atmospheric_forecast_runs_cold(predictions):
    """Directional prediction tied to the absent ENSO term."""
    atm = trend_comparison(predictions)['atmospheric']
    assert atm['forecast_per_year'] < atm['observed_per_year']


# --- the record ------------------------------------------------------------

def test_falsification_criteria_cover_every_system(predictions):
    text = " ".join(falsification_criteria(predictions))
    for name in SYSTEMS:
        assert name in text
    assert "whole model" in text


def test_record_is_serialisable_and_complete(predictions):
    import json
    rec = to_record(predictions)
    json.dumps(rec)  # must not raise

    assert rec["initialised_from"] == LAST_OBSERVED_YEAR
    assert rec["model_has_enso"] is False
    assert len(rec["predictions"]) == len(predictions)
    assert rec["falsification"]
    assert set(rec["resolution_sources"]) == set(SYSTEMS)


def test_record_states_the_interval_basis(predictions):
    """An interval with no stated basis is not pre-registration."""
    rec = to_record(predictions)
    assert "held-out" in rec["interval_basis"]


def test_forecast_is_reproducible(priors):
    """Same inputs, same numbers — a forecast that drifts cannot be checked."""
    a = run_forecast(priors=priors)
    b = run_forecast(priors=priors)
    assert [p.value for p in a] == pytest.approx([p.value for p in b])
