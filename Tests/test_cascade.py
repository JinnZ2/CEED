# Tests/test_cascade.py

"""
Tests for coupled tipping elements.

The load-bearing ones are `test_uncoupled_trigger_does_not_spread` and
`test_coupling_propagates_the_tip` — together they show the cascade is caused
by the coupling and not by the forcing. `test_influence_arrives_through_the_
slow_variable` pins the design choice that makes cascades delay-dependent.
"""

import numpy as np
import pytest

from simulation.cascade import (
    NetworkElement,
    TippingNetwork,
    thwaites_ross_network,
)
from simulation.tipping import FOLD_F


# Shorter horizon / coarser step for the bisection-heavy threshold tests.
FAST = dict(years=2000.0, dt=0.5)


def two_node(c_ab: float = 0.0, c_ba: float = 0.0) -> TippingNetwork:
    """Minimal A -> B network with symmetric timescales."""
    return TippingNetwork(
        elements=[
            NetworkElement("a", tau_fast=1.0, tau_slow=40.0, x0=-1.0),
            NetworkElement("b", tau_fast=1.0, tau_slow=40.0, x0=-1.0),
        ],
        coupling=np.array([[0.0, c_ba], [c_ab, 0.0]]),
    )


# --- construction ----------------------------------------------------------

def test_rejects_wrong_shaped_coupling():
    with pytest.raises(ValueError):
        TippingNetwork(elements=[NetworkElement("a")],
                       coupling=np.zeros((2, 2)))


def test_rejects_nonzero_diagonal():
    with pytest.raises(ValueError):
        two_node().__class__(
            elements=[NetworkElement("a"), NetworkElement("b")],
            coupling=np.array([[0.5, 0.0], [0.0, 0.0]]))


def test_element_rejects_inverted_timescales():
    with pytest.raises(ValueError):
        NetworkElement("bad", tau_fast=10.0, tau_slow=1.0)


def test_state_starts_at_initial_conditions():
    net = two_node()
    x, h = net.split()
    assert np.allclose(x, [-1.0, -1.0])
    assert np.allclose(h, [-1.0, -1.0])


# --- the cascade itself ----------------------------------------------------

def test_uncoupled_trigger_does_not_spread():
    """With no coupling, pushing A must leave B exactly where it was."""
    net = two_node(c_ab=0.0)
    times = net.cascade_from("a")
    assert times["a"] is not None
    assert times["b"] is None


def test_coupling_propagates_the_tip():
    """Same push, coupling on: B tips although it is never forced."""
    net = two_node(c_ab=0.8)
    times = net.cascade_from("a")
    assert times["a"] is not None
    assert times["b"] is not None
    assert times["b"] > times["a"]


def test_cascade_survives_removal_of_the_trigger():
    """
    The push lasts 20 units. If B tips after that, the cascade is being
    driven by A's state, not by any residual forcing.
    """
    net = two_node(c_ab=0.8)
    times = net.cascade_from("a", push_years=20.0, years=4000.0)
    assert times["b"] > 20.0


def test_stabilising_coupling_prevents_the_cascade():
    """Negative coupling is protective; not every interaction destabilises."""
    net = two_node(c_ab=-0.8)
    assert net.cascade_from("a")["b"] is None


def test_cascade_from_rejects_unknown_element():
    with pytest.raises(ValueError):
        two_node().cascade_from("nonexistent")


def test_domino_threshold_is_bracketed_and_reproducible():
    net = two_node()
    first = net.domino_threshold("a", "b", **FAST)
    second = net.domino_threshold("a", "b", **FAST)

    assert first is not None
    assert 0.0 < first < 2.0
    assert first == pytest.approx(second)


def test_domino_threshold_restores_the_matrix():
    net = two_node(c_ab=0.25)
    net.domino_threshold("a", "b", **FAST)
    assert net.coupling[1, 0] == pytest.approx(0.25)


def test_below_threshold_no_cascade_above_threshold_cascade():
    net = two_node()
    thr = net.domino_threshold("a", "b", **FAST)

    below = two_node(c_ab=thr * 0.9)
    above = two_node(c_ab=thr * 1.1)
    assert below.cascade_from("a")["b"] is None
    assert above.cascade_from("a")["b"] is not None


def test_no_cascade_returns_none_threshold():
    """A purely stabilising link has no domino threshold."""
    net = two_node()
    assert net.domino_threshold("a", "b", lo=-2.0, hi=-0.1, **FAST) is None


# --- the design choice being tested ----------------------------------------

def test_influence_arrives_through_the_slow_variable():
    """
    Coupling reads h, not x. Two states with identical x but different h must
    produce different derivatives — that is what makes cascades delayed.
    """
    net = two_node(c_ab=0.8)
    same_x = -1.0

    early = np.array([same_x, same_x, -1.0, -1.0])   # h not yet responded
    late = np.array([same_x, same_x, +1.0, -1.0])    # h has responded

    d_early = net.derivative(early, np.zeros(2))
    d_late = net.derivative(late, np.zeros(2))

    assert d_early[1] != pytest.approx(d_late[1])


def test_untouched_network_exerts_no_influence():
    """
    Regression. Coupling reads the tipped fraction (h+1)/2, not raw h. Reading
    h directly gave every element a shove at t=0 (since h starts at -1) and
    inverted stabilising links, so a negative coefficient tipped its target.
    """
    net = two_node(c_ab=-0.8, c_ba=-0.8)
    intact = np.array([-1.0, -1.0, -1.0, -1.0])
    d = net.derivative(intact, np.zeros(2))

    # x = -1 is an equilibrium of x - x^3, so with zero influence dx must be 0.
    assert d[0] == pytest.approx(0.0, abs=1e-12)
    assert d[1] == pytest.approx(0.0, abs=1e-12)


def test_stabilising_link_is_protective_not_inverted():
    """A negative coefficient must push the target away from its fold."""
    net = two_node(c_ab=-0.8)
    tipped_source = np.array([-1.0, -1.0, +1.0, -1.0])   # a fully tipped
    d = net.derivative(tipped_source, np.zeros(2))
    assert d[1] < 0          # b is pushed further down, not up


def test_slower_response_delays_the_cascade():
    """Larger tau_slow on the trigger means its influence arrives later."""
    delays = []
    for tau_slow in (20.0, 200.0):
        net = TippingNetwork(
            elements=[
                NetworkElement("a", tau_fast=1.0, tau_slow=tau_slow, x0=-1.0),
                NetworkElement("b", tau_fast=1.0, tau_slow=40.0, x0=-1.0),
            ],
            coupling=np.array([[0.0, 0.0], [0.9, 0.0]]),
        )
        times = net.cascade_from("a", years=6000.0)
        delays.append(times["b"] - times["a"])
    assert delays[1] > delays[0]


def test_asymmetric_coupling_is_directional():
    """A -> B present, B -> A absent: tipping B must not tip A."""
    net = two_node(c_ab=0.8, c_ba=0.0)
    assert net.cascade_from("b")["a"] is None


# --- diagnostics -----------------------------------------------------------

def test_basin_margin_orders_elements_by_proximity():
    net = TippingNetwork(
        elements=[NetworkElement("near", x0=-0.2),
                  NetworkElement("far", x0=-1.0)],
        coupling=np.zeros((2, 2)),
    )
    margins = net.basin_margin()
    assert margins[0] < margins[1]
    assert margins[0] == pytest.approx(0.2)
    assert margins[1] == pytest.approx(1.0)


def test_basin_margin_is_negative_once_past_the_boundary():
    net = two_node()
    assert net.basin_margin(np.array([0.5, -1.0]))[0] < 0


def test_basin_margin_beyond_the_fold_has_no_boundary():
    net = two_node()
    m = net.basin_margin(np.array([-1.0, -1.0]),
                         forcings=np.array([FOLD_F * 1.5, 0.0]))
    assert m[0] == -np.inf


def test_consequence_is_clipped_at_full():
    """Destabilising coupling pushes h past +1; that is still 100%, not more."""
    net = TippingNetwork(
        elements=[NetworkElement("a", consequence=2.0),
                  NetworkElement("b", consequence=3.0)],
        coupling=np.zeros((2, 2)),
    )
    assert net.consequence_realised(np.array([5.0, 5.0])) == pytest.approx(5.0)
    assert net.consequence_realised(np.array([-1.0, -1.0])) == pytest.approx(0.0)
    assert net.consequence_realised(np.array([1.0, -1.0])) == pytest.approx(2.0)


def test_divergence_fails_loudly():
    """A too-large step must raise, not silently emit NaN. 'Never tipped' and
    'blew up' must not look alike in the diagnostics."""
    net = two_node(c_ab=0.5)
    # Overflow is expected on the way to the guard; it is the point of the test.
    with np.errstate(over='ignore', invalid='ignore'):
        with pytest.raises(FloatingPointError):
            net.simulate(lambda tt: np.array([0.45, 0.0]), years=200, dt=3.0)


def test_tipped_mask():
    net = two_node()
    assert list(net.tipped(np.array([0.5, -0.5]))) == [True, False]


# --- determinism and shape -------------------------------------------------

def test_derivative_is_deterministic():
    net = two_node(c_ab=0.5)
    s = np.array([-0.3, 0.2, -0.5, 0.1])
    F = np.array([0.1, 0.2])
    nz = np.array([0.01, -0.01])
    assert np.array_equal(net.derivative(s, F, nz), net.derivative(s, F, nz))


def test_simulate_rejects_short_noise():
    net = two_node()
    with pytest.raises(ValueError):
        net.simulate(lambda t: np.zeros(2), years=100, dt=0.5,
                     noise=np.zeros((5, 2)))


def test_trajectory_shape_and_finiteness():
    net = two_node(c_ab=0.5)
    t, x, h = net.simulate(lambda tt: np.zeros(2), years=100, dt=0.5)
    assert x.shape == h.shape == (len(t), 2)
    assert np.all(np.isfinite(x)) and np.all(np.isfinite(h))


# --- the illustrative configuration ----------------------------------------

def test_thwaites_ross_encodes_the_asymmetry():
    """
    The element nearer its threshold carries the smaller consequence. That
    inversion is the point of the configuration.
    """
    net = thwaites_ross_network()
    margins = net.basin_margin()
    thw, ross = net.elements

    assert margins[0] < margins[1]              # thwaites is closer
    assert thw.consequence < ross.consequence   # and matters far less
    assert ross.consequence / thw.consequence > 10


def test_thwaites_ross_cascade_is_caused_by_coupling():
    uncoupled = thwaites_ross_network(coupling_strength=0.0)
    coupled = thwaites_ross_network(coupling_strength=0.55)

    assert uncoupled.cascade_from("thwaites")["ross"] is None
    assert coupled.cascade_from("thwaites")["ross"] is not None


def test_most_consequence_arrives_after_the_trigger_is_gone():
    """Thwaites commits early and cheaply; Ross carries the cost, late."""
    net = thwaites_ross_network(coupling_strength=0.55)
    t, x, h = net.simulate(
        lambda tt: np.array([0.45 if tt < 100 else 0.0, 0.0]),
        years=4000, dt=0.5)

    total = sum(e.consequence for e in net.elements)
    at_trigger_end = net.consequence_realised(h[int(100 / 0.5)])
    at_end = net.consequence_realised(h[-1])

    assert at_trigger_end / total < 0.15
    assert at_end / total > 0.95
