# Tests/test_validation.py

"""
Tests for the held-out validation split.

The load-bearing test here is `test_calibration_cannot_see_the_test_window`:
it corrupts the test-window observations and asserts the fitted parameters do
not move. Everything else about a train/test split is decoration if that
property does not hold.
"""

import numpy as np
import pytest

from simulation.convergence_model import AnthropogenicForcing, SystemParameters
from simulation.unit_bridge import ALL_OBSERVED
from simulation.validation import (
    CALIBRATION_TARGETS,
    DEFAULT_END,
    DEFAULT_START,
    DEFAULT_TRAIN_END,
    calibrate,
    derive_solar_modulation,
    evaluate_split,
)

# Small search settings so the suite stays fast.
FAST = dict(refine_rounds=1, grid=4)


# --- derived constants -----------------------------------------------------

def test_solar_modulation_is_derived_not_hardcoded():
    params = SystemParameters()
    value = derive_solar_modulation(params, 2010, 2024)
    # Matches the shipped constant, which was derived the same way.
    assert value == pytest.approx(params.solar_modulation, rel=0.02)


def test_train_window_gives_a_different_modulation():
    """
    The full-record F10.7 range includes the 2024 maximum of 180 sfu, which
    sits in the test window. A train-only derivation must therefore differ —
    if it did not, the constant would be carrying test information.
    """
    params = SystemParameters()
    full = derive_solar_modulation(params, DEFAULT_START, DEFAULT_END)
    train = derive_solar_modulation(params, DEFAULT_START, DEFAULT_TRAIN_END)

    assert train < full
    assert train == pytest.approx(0.337, abs=0.01)


def test_modulation_scales_with_observed_amplitude():
    params = SystemParameters()
    wide = derive_solar_modulation(params, 2010, 2024)   # range 69-180
    narrow = derive_solar_modulation(params, 2017, 2020)  # range 69-77
    assert narrow < wide


def test_modulation_rejects_a_degenerate_window():
    with pytest.raises(ValueError):
        derive_solar_modulation(SystemParameters(), 2015, 2015)


# --- the leakage guarantee -------------------------------------------------

def test_calibration_cannot_see_the_test_window(monkeypatch):
    """
    Corrupt every test-window observation and refit. If the calibration is
    genuinely confined to the training window, the parameters are unchanged.
    """
    before = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)

    for system, table in ALL_OBSERVED.items():
        poisoned = dict(table)
        for year in range(DEFAULT_TRAIN_END + 1, DEFAULT_END + 1):
            if year in poisoned:
                poisoned[year] = poisoned[year] * 10.0 + 1234.0
        monkeypatch.setitem(ALL_OBSERVED, system, poisoned)

    after = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)

    assert after.A_0 == pytest.approx(before.A_0)
    assert after.atm_fraction == pytest.approx(before.atm_fraction)
    assert after.solar_modulation == pytest.approx(before.solar_modulation)


def test_calibration_does_respond_to_training_data(monkeypatch):
    """
    The mirror of the previous test: poisoning the TRAIN window must move the
    fit. Otherwise the leakage test above would pass trivially.
    """
    before = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)

    table = dict(ALL_OBSERVED['oceanic'])
    for year in range(DEFAULT_START, DEFAULT_TRAIN_END + 1):
        if year in table:
            table[year] = table[year] * 4.0
    monkeypatch.setitem(ALL_OBSERVED, 'oceanic', table)

    after = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)

    assert (after.A_0 != pytest.approx(before.A_0)
            or after.atm_fraction != pytest.approx(before.atm_fraction))


# --- calibration behaviour -------------------------------------------------

def test_calibration_is_deterministic():
    """Grid search, not a stochastic optimiser — same window, same answer."""
    first = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)
    second = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)

    assert first.A_0 == pytest.approx(second.A_0)
    assert first.atm_fraction == pytest.approx(second.atm_fraction)
    assert first.train_objective == pytest.approx(second.train_objective)


def test_calibration_returns_usable_parameters():
    result = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)

    assert isinstance(result.params, SystemParameters)
    assert result.params.anthropogenic is not None
    assert result.params.anthropogenic.A_0 == pytest.approx(result.A_0)
    assert result.params.anthropogenic.atm_fraction == pytest.approx(
        result.atm_fraction)
    assert np.isfinite(result.train_objective)
    assert result.n_evaluations > 0


def test_calibration_beats_an_arbitrary_parameter_choice():
    """A fit that cannot beat a deliberately poor guess is not fitting."""
    from simulation.validation import _objective

    result = calibrate(train_end=DEFAULT_TRAIN_END, **FAST)
    poor = SystemParameters(
        anthropogenic=AnthropogenicForcing(A_0=6.0, atm_fraction=0.6))

    assert result.train_objective < _objective(
        poor, DEFAULT_START, DEFAULT_TRAIN_END)


# --- split evaluation ------------------------------------------------------

def test_split_windows_are_disjoint_and_complete():
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    split = evaluate_split(params, "test", systems=list(CALIBRATION_TARGETS))

    train_years = set(split.train['oceanic'].years)
    test_years = set(split.test['oceanic'].years)

    assert train_years.isdisjoint(test_years)
    assert max(train_years) == DEFAULT_TRAIN_END
    assert min(test_years) == DEFAULT_TRAIN_END + 1
    assert max(test_years) == DEFAULT_END


def test_initialisation_year_is_in_neither_window():
    """2010 is fit by construction, so it must not be scored at all."""
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    split = evaluate_split(params, "test", systems=['oceanic'])

    assert DEFAULT_START not in split.train['oceanic'].years
    assert DEFAULT_START not in split.test['oceanic'].years


def test_test_window_is_a_single_trajectory_not_a_restart():
    """
    The test window must be a free-running continuation, not a warm start.
    Re-initialising at 2020 would make the first test year near-perfect; a
    genuine 10-year-old forecast will not be.
    """
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    split = evaluate_split(params, "test", systems=['oceanic'])

    first_test_error = abs(split.test['oceanic'].predicted[0]
                           - split.test['oceanic'].observed[0])
    assert first_test_error > 1e-6


def test_gap_is_train_minus_test():
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    split = evaluate_split(params, "test", systems=['oceanic'])

    expected = (split.train['oceanic'].r_squared
                - split.test['oceanic'].r_squared)
    assert split.gap('oceanic') == pytest.approx(expected)


def test_mean_r2_matches_the_scored_systems():
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    systems = list(CALIBRATION_TARGETS)
    split = evaluate_split(params, "test", systems=systems)

    expected = np.mean([split.test[s].r_squared for s in systems])
    assert split.mean_r2('test') == pytest.approx(expected)


def test_shipped_parameters_overfit_the_ocean():
    """
    Records the headline result: oceanic scores well in-sample and fails
    out-of-sample under the shipped parameters. If this ever stops being
    true the audit and CLAUDE.md need updating with it.
    """
    params = SystemParameters(anthropogenic=AnthropogenicForcing())
    split = evaluate_split(params, "test", systems=['oceanic'])

    assert split.train['oceanic'].r_squared > 0.5
    assert split.test['oceanic'].r_squared < 0.0
    assert split.gap('oceanic') > 1.0
