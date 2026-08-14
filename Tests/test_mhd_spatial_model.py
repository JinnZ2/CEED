# Tests/test_mhd_spatial_model.py

import numpy as np
import pytest

from simulation.mhd_spatial_model import (
    COUPLING_MATRIX,
    EARTH_OMEGA,
    MU_0,
    ZONE_AREA_FRACTION,
    ZONES,
    IonosphereBaseline,
    SolarWindBaseline,
    SpatialCEED,
    Zone,
    current_density,
    dynamo_sheet_current,
    effective_torque,
    ionospheric_dynamo_power_density,
    ionospheric_dynamo_spec,
    lorentz_force_density,
    magnetic_torque,
    mhd_energy_injection,
    motional_electric_field,
    poynting_flux,
    solar_wind_kinetic_flux,
)


# --- constants -------------------------------------------------------------

def test_physical_constants():
    assert MU_0 == pytest.approx(1.25663706e-6, rel=1e-6)
    assert EARTH_OMEGA == pytest.approx(2 * np.pi / 86164.0, rel=1e-4)


def test_solar_wind_mass_density():
    sw = SolarWindBaseline()
    # 5 protons/cm^3 -> ~8.4e-21 kg/m^3
    assert sw.mass_density == pytest.approx(8.36e-21, rel=1e-2)


# --- MHD injection ---------------------------------------------------------

def test_motional_field_is_perpendicular():
    v = [4.0e5, 0.0, 0.0]
    B = [0.0, 0.0, 5.0e-9]
    E = motional_electric_field(v, B)

    assert np.dot(E, v) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(E, B) == pytest.approx(0.0, abs=1e-12)
    assert np.linalg.norm(E) == pytest.approx(4.0e5 * 5.0e-9)


def test_parallel_flow_injects_nothing():
    """v parallel to B gives v x B = 0, so no injection."""
    v = [1.0e5, 0.0, 0.0]
    B = [3.0e-9, 0.0, 0.0]
    assert mhd_energy_injection(v, B, rho=1e-20) == pytest.approx(0.0)


def test_injection_scales_quadratically_with_speed():
    B = [0.0, 0.0, 5.0e-9]
    rho = 8.36e-21
    single = mhd_energy_injection([1.0e5, 0, 0], B, rho)
    double = mhd_energy_injection([2.0e5, 0, 0], B, rho)
    assert double == pytest.approx(4.0 * single)


def test_injection_inversely_proportional_to_density():
    v = [4.0e5, 0.0, 0.0]
    B = [0.0, 0.0, 5.0e-9]
    assert (mhd_energy_injection(v, B, rho=2e-20)
            == pytest.approx(0.5 * mhd_energy_injection(v, B, rho=1e-20)))


def test_injection_rejects_nonpositive_density():
    with pytest.raises(ValueError):
        mhd_energy_injection([1.0, 0, 0], [0, 0, 1.0], rho=0.0)


def test_solar_wind_energy_fluxes_are_physical():
    """
    Nominal solar wind at 1 AU carries a bulk kinetic energy flux of order
    1e-4 W/m^2, and its Poynting flux is one to two orders smaller.
    """
    sw = SolarWindBaseline()
    kinetic = solar_wind_kinetic_flux(sw.speed, sw.mass_density)
    assert 1e-5 < kinetic < 1e-3

    S = poynting_flux([sw.speed, 0, 0], [0, 0, sw.imf_magnitude])
    assert 0 < S < kinetic


# --- torque ----------------------------------------------------------------

def test_uniform_field_carries_no_current():
    """A spatially uniform B has zero curl, hence zero current."""
    grid = np.zeros((3, 5, 5, 5))
    grid[2] = 4.5e-5  # uniform Bz
    J = current_density(grid, spacing=1.0e3)
    assert np.allclose(J, 0.0, atol=1e-18)


def test_current_density_recovers_known_curl():
    """
    For B = (0, b*x, 0), curl(B) = (0, 0, b), so J_z = b / mu_0.
    """
    n, spacing, b = 6, 2.0, 3.0e-9
    x = np.arange(n) * spacing
    grid = np.zeros((3, n, n, n))
    grid[1] = x[:, None, None]* b

    J = current_density(grid, spacing=spacing)
    interior = J[2, 1:-1, 1:-1, 1:-1]
    assert np.allclose(interior, b / MU_0, rtol=1e-9)


def test_current_density_rejects_bad_shape():
    with pytest.raises(ValueError):
        current_density(np.zeros((2, 4, 4, 4)), spacing=1.0)


def test_lorentz_force_perpendicular_to_both():
    J = np.array([1.0, 0.0, 0.0])
    B = np.array([0.0, 2.0, 0.0])
    F = lorentz_force_density(J, B)

    assert np.dot(F, J) == pytest.approx(0.0)
    assert np.dot(F, B) == pytest.approx(0.0)
    assert F[2] == pytest.approx(2.0)


def test_force_parallel_to_field_vanishes():
    J = np.array([0.0, 0.0, 3.0])
    B = np.array([0.0, 0.0, 7.0])
    assert np.allclose(lorentz_force_density(J, B), 0.0)


def test_magnetic_torque_of_radial_force_is_zero():
    """r x F vanishes when F is parallel to r."""
    r = np.array([[1.0], [0.0], [0.0]])
    J = np.array([[1.0], [0.0], [0.0]])
    B = np.array([[2.0], [0.0], [0.0]])  # J x B = 0
    assert np.allclose(magnetic_torque(r, J, B, cell_volume=1.0), 0.0)


def test_magnetic_torque_known_value():
    """r = x_hat, J = y_hat, B = z_hat -> J x B = x_hat, r x F = 0."""
    r = np.array([[0.0], [1.0], [0.0]])
    J = np.array([[0.0], [1.0], [0.0]])
    B = np.array([[0.0], [0.0], [1.0]])
    # J x B = x_hat ; r x F = y_hat x x_hat = -z_hat
    tau = magnetic_torque(r, J, B, cell_volume=2.0)
    assert tau == pytest.approx([0.0, 0.0, -2.0])


def test_effective_torque_is_antisymmetric_in_latitude():
    kwargs = dict(lon_deg=0.0, B_magnitude=4.5e-5, v_magnitude=250.0)
    north = effective_torque(lat_deg=45.0, **kwargs)
    south = effective_torque(lat_deg=-45.0, **kwargs)

    assert north == pytest.approx(-south)
    assert north > 0


def test_effective_torque_vanishes_at_equator():
    tau = effective_torque(0.0, 0.0, 4.5e-5, 250.0)
    assert tau == pytest.approx(0.0)


# --- dynamo ----------------------------------------------------------------

def test_dynamo_power_is_positive_and_quadratic():
    v = [250.0, 0.0, 0.0]
    B = [0.0, 0.0, 4.5e-5]
    single = ionospheric_dynamo_power_density(v, B)
    double = ionospheric_dynamo_power_density([500.0, 0.0, 0.0], B)

    assert single > 0
    assert double == pytest.approx(4.0 * single)


def test_dynamo_joule_heating_magnitude_is_physical():
    """
    Height-integrated Joule dissipation in the dynamo region should land in
    the milliwatt-per-square-metre band, not kilowatts.
    """
    iono = IonosphereBaseline()
    P = ionospheric_dynamo_power_density(
        [iono.neutral_wind_speed, 0.0, 0.0],
        [0.0, 0.0, iono.field_magnitude],
        layer_height=iono.layer_height,
        conductivity=iono.conductivity,
    )
    assert 1e-4 < P < 1.0


def test_dynamo_sheet_current_relates_to_power():
    """P = K * |v x B| for the same inputs."""
    v = [250.0, 0.0, 0.0]
    B = [0.0, 0.0, 4.5e-5]
    K = dynamo_sheet_current(v, B)
    E = np.linalg.norm(motional_electric_field(v, B))
    assert ionospheric_dynamo_power_density(v, B) == pytest.approx(K * E)


def test_dynamo_spec_formula_matches_its_definition():
    v = [250.0, 0.0, 0.0]
    B = [0.0, 0.0, 4.5e-5]
    expected = EARTH_OMEGA * dynamo_sheet_current(v, B)
    assert ionospheric_dynamo_spec(v, B) == pytest.approx(expected)


# --- spatial zoning --------------------------------------------------------

def test_latitude_bands_tile_the_sphere():
    bands = sum(ZONE_AREA_FRACTION[z]
                for z in (Zone.POLAR, Zone.MIDLAT, Zone.EQUATOR))
    assert bands == pytest.approx(1.0)


def test_coupling_matrix_is_symmetric_with_zero_diagonal():
    assert np.allclose(COUPLING_MATRIX, COUPLING_MATRIX.T)
    assert np.allclose(np.diag(COUPLING_MATRIX), 0.0)


def test_exchange_conserves_energy():
    model = SpatialCEED(energy=np.array([300.0, 50.0, 100.0, 20.0]))
    assert model.exchange(model.energy).sum() == pytest.approx(0.0, abs=1e-12)


def test_closed_system_conserves_total_energy():
    """No sources and no sinks: coupling only redistributes."""
    model = SpatialCEED(energy=np.array([300.0, 50.0, 100.0, 20.0]))
    start = model.energy.sum()
    _, history = model.run(years=10.0)
    assert history[-1].sum() == pytest.approx(start, rel=1e-10)


def test_coupling_drives_zones_toward_equality():
    model = SpatialCEED(energy=np.array([300.0, 50.0, 100.0, 20.0]))
    spread_before = model.energy.max() - model.energy.min()
    _, history = model.run(years=20.0)
    spread_after = history[-1].max() - history[-1].min()
    assert spread_after < spread_before


def test_energy_flows_from_high_to_low():
    model = SpatialCEED(energy=np.array([500.0, 0.0, 0.0, 0.0]))
    flow = model.exchange(model.energy)
    assert flow[0] < 0           # the hot zone loses
    assert all(f > 0 for f in flow[1:])  # the others gain


def test_sinks_remove_energy():
    model = SpatialCEED(
        energy=np.array([100.0, 100.0, 100.0, 100.0]),
        sink_rates=np.full(4, 0.1),
    )
    _, history = model.run(years=5.0)
    assert history[-1].sum() < history[0].sum()


def test_history_shape_and_time_axis():
    model = SpatialCEED()
    t, history = model.run(years=3.0, dt=1.0 / 12.0)
    assert history.shape == (len(t), len(ZONES))
    assert t[0] == 0.0
    assert t[-1] == pytest.approx(3.0)


def test_as_dict_keys_match_zones():
    model = SpatialCEED()
    d = model.as_dict()
    assert set(d) == {z.value for z in ZONES}


def test_rejects_nonzero_diagonal_coupling():
    bad = COUPLING_MATRIX.copy()
    bad[0, 0] = 0.5
    with pytest.raises(ValueError):
        SpatialCEED(coupling=bad)


def test_rejects_wrong_energy_length():
    with pytest.raises(ValueError):
        SpatialCEED(energy=np.array([1.0, 2.0]))
