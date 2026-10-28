"""Tests for simulation/radiation_belt.py (barrier slot, event-induced crossing)."""
import ast
import inspect

import numpy as np
import pytest

from simulation import radiation_belt as rb
from simulation.convergence_model import ConvergencePredictor
from simulation.tipping import ForcingEvent

T = np.arange(0.0, 400.0, 1.0)
QUIET = np.full_like(T, 2.0)


def test_no_crossing_below_barrier_inner_stays_zero():
    r = rb.run_belt(QUIET, T, [ForcingEvent(50.0, 2.0, 2.0)], barrier_kp=6.0)
    assert not r.crossed
    assert np.all(r.inner == 0.0)
    assert r.outer.max() > 0.0  # the outer belt still fills


def test_one_storm_above_barrier_crosses_once_and_fills_inner():
    r = rb.run_belt(QUIET, T, [ForcingEvent(50.0, 5.0, 2.0)], barrier_kp=6.0)
    assert r.crossing_days == [50.0]
    assert r.inner.max() > 0.0


def test_inner_decays_with_carried_efolding():
    # start with inner filled, nothing driving it
    r = rb.run_belt(np.zeros_like(T), T, [], barrier_kp=6.0, inner0=1.0)
    k = int(rb.INNER_EFOLD_DAYS)
    expected = (1 - 1 / rb.INNER_EFOLD_DAYS) ** k  # explicit Euler, dt = 1 d
    assert r.inner[k] == pytest.approx(expected, rel=1e-12)
    assert r.inner[k] == pytest.approx(np.exp(-1.0), abs=0.01)


def test_barrier_is_a_reading_of_its_window():
    window = [2.0, 4.5, 6.2, 3.0]
    b = rb.characterized_barrier(window)
    assert b == 6.2
    # mean driver far below the barrier, one later storm exceeds it
    r = rb.run_belt(QUIET, T, [ForcingEvent(100.0, 4.5, 1.0)], barrier_kp=b)
    assert QUIET.mean() < b and r.crossed


def test_kp_ceiling_saturates_and_barrier_at_ceiling_cannot_be_crossed():
    b = rb.characterized_barrier([2.0, 12.0])  # an over-scale storm in the window
    assert b == rb.KP_MAX and rb.barrier_reading(b) == "AT_CEILING"
    r = rb.run_belt(QUIET, T, [ForcingEvent(100.0, 30.0, 3.0)], barrier_kp=b)
    assert r.kp.max() == rb.KP_MAX
    assert r.saturated_days == 3
    assert not r.crossed  # the driver has no reading above its own ceiling
    assert rb.barrier_reading(6.0) == "BELOW_CEILING"


def test_empty_window_rejected():
    with pytest.raises(ValueError):
        rb.characterized_barrier([])


def test_storm_train_deterministic_and_in_days():
    a = rb.storm_train(4.0, 3.0, 1.5, seed=3)
    b = rb.storm_train(4.0, 3.0, 1.5, seed=3)
    assert [(e.time, e.magnitude) for e in a] == [(e.time, e.magnitude) for e in b]
    assert all(0.0 <= e.time < 3.0 * rb.DAYS_PER_YEAR for e in a)
    assert all(e.magnitude >= 0.0 for e in a)


def test_run_belt_draws_no_random_numbers():
    tree = ast.parse(inspect.getsource(rb.run_belt))
    names = {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    assert not names & {"random", "normal", "uniform", "poisson", "exponential",
                        "default_rng"}


def test_wired_into_predictor_without_changing_default_output():
    p = ConvergencePredictor()
    t, sol = p.predict_convergence(years=3.0)
    t2, sol2 = ConvergencePredictor().predict_convergence(years=3.0)
    assert np.array_equal(sol, sol2)
    r = p.radiation_belt(t, sol, rb.storm_train(4.0, 3.0, 1.5, seed=0), 6.0)
    assert r.t_days[-1] == pytest.approx(t[-1] * rb.DAYS_PER_YEAR, abs=1.0)
    assert np.all(np.isfinite(r.inner)) and np.all(r.kp <= rb.KP_MAX)
    assert np.array_equal(sol, sol2)  # reading the belt changed nothing


def test_mismatched_inputs_rejected():
    with pytest.raises(ValueError):
        rb.run_belt(QUIET[:-1], T, [], 6.0)
