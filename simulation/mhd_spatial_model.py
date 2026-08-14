# simulation/mhd_spatial_model.py

"""
CEED MHD and Spatial Zoning Module

Implements the four roadmap items from To-be-added.md:

  1. MHD energy injection      E_input_mhd(t) = eta * |v_plasma x B|^2 / rho
  2. Magnetic feedback torque  tau = integral( r x (J x B) ) dV
                               tau_eff = beta * sin(lat) * cos(lon) * |B| * |v|
  3. Spatial zoning            4 zones with a coupling matrix, replacing the
                               single-point ("sad little dot") Earth
  4. Earth-ionosphere dynamo   E_dynamo(t) = omega_earth * (v_atmo x B) * h * sigma

UNITS NOTICE
------------
Two of the roadmap formulas are not dimensionally closed as written. They are
implemented here verbatim so the specification is honoured, and each is paired
with a dimensionally consistent companion whose units are stated explicitly:

  - `mhd_energy_injection` (spec) vs `poynting_flux` / `solar_wind_kinetic_flux`
  - `ionospheric_dynamo_spec` (spec) vs `ionospheric_dynamo_power_density`

Prefer the companions for anything quantitative. See Docs/numerical-audit.md.

All physical constants are SI unless a name says otherwise.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

MU_0 = 4.0e-7 * np.pi              # vacuum permeability, H/m (T*m/A)
PROTON_MASS = 1.67262192369e-27    # kg
EARTH_RADIUS = 6.371e6             # m (mean)
EARTH_OMEGA = 7.2921150e-5         # rad/s (sidereal rotation rate)


@dataclass(frozen=True)
class SolarWindBaseline:
    """
    Nominal quiet-time solar wind at 1 AU.

    Values are the standard nominal set used in magnetospheric modelling:
    ~400 km/s bulk speed, ~5 protons/cm^3, ~5 nT interplanetary field.
    Storm-time values differ by roughly an order of magnitude in density
    and a factor of 2-3 in speed.
    """
    speed: float = 4.0e5            # m/s
    number_density: float = 5.0e6   # m^-3  (5 cm^-3)
    imf_magnitude: float = 5.0e-9   # T     (5 nT)

    @property
    def mass_density(self) -> float:
        """Proton mass density, kg/m^3."""
        return self.number_density * PROTON_MASS


@dataclass(frozen=True)
class IonosphereBaseline:
    """
    Nominal Earth-ionosphere dynamo region parameters.

    Note on `layer_height`: the roadmap specifies 110 km. Hall conductivity
    peaks near 110 km while Pedersen conductivity peaks nearer 125 km, so
    110 km is the Hall-weighted choice. `conductivity` spans the range given
    in the roadmap (1e-4 to 1e-2 S/m); the default sits at the geometric
    middle of that span.
    """
    layer_height: float = 1.10e5        # m (110 km)
    conductivity: float = 1.0e-3        # S/m
    conductivity_range: Tuple[float, float] = (1.0e-4, 1.0e-2)  # S/m
    neutral_wind_speed: float = 250.0   # m/s (roadmap value; storm-time upper end)
    field_magnitude: float = 4.5e-5     # T (45,000 nT, mid-range surface field)
    field_range: Tuple[float, float] = (2.5e-5, 6.5e-5)  # T (25,000-65,000 nT)


# ---------------------------------------------------------------------------
# 1. MHD energy injection
# ---------------------------------------------------------------------------

def motional_electric_field(v_plasma: Sequence[float],
                            B: Sequence[float]) -> np.ndarray:
    """
    Motional electric field E = v x B.

    Args:
        v_plasma: plasma velocity vector, m/s
        B: magnetic field vector, T

    Returns:
        Electric field vector, V/m
    """
    return np.cross(np.asarray(v_plasma, dtype=float),
                    np.asarray(B, dtype=float))


def mhd_energy_injection(v_plasma: Sequence[float],
                         B: Sequence[float],
                         rho: float,
                         eta: float = 1.0e-3) -> float:
    """
    Roadmap formula, implemented verbatim:

        E_input_mhd = eta * |v_plasma x B|^2 / rho

    Args:
        v_plasma: plasma velocity vector, m/s
        B: magnetic field vector, T
        rho: plasma mass density, kg/m^3 (must be > 0)
        eta: plasma conductivity scaling, roadmap range ~0.001-1

    Returns:
        Injection term in the roadmap's own (non-SI-closed) units.
        See `poynting_flux` for a W/m^2 quantity.

    Raises:
        ValueError: if rho is not strictly positive.
    """
    if rho <= 0:
        raise ValueError(f"plasma mass density must be > 0, got {rho}")

    E = motional_electric_field(v_plasma, B)
    return eta * float(np.dot(E, E)) / rho


def poynting_flux(v_plasma: Sequence[float], B: Sequence[float]) -> float:
    """
    Dimensionally closed companion to `mhd_energy_injection`.

    For a perfectly conducting flow, S = |(v x B) x B| / mu_0.

    Args:
        v_plasma: plasma velocity vector, m/s
        B: magnetic field vector, T

    Returns:
        Electromagnetic energy flux, W/m^2
    """
    E = motional_electric_field(v_plasma, B)
    S = np.cross(E, np.asarray(B, dtype=float)) / MU_0
    return float(np.linalg.norm(S))


def solar_wind_kinetic_flux(speed: float, mass_density: float) -> float:
    """
    Bulk kinetic energy flux of the solar wind, 0.5 * rho * v^3.

    Args:
        speed: bulk speed, m/s
        mass_density: kg/m^3

    Returns:
        Kinetic energy flux, W/m^2
    """
    return 0.5 * mass_density * speed**3


# ---------------------------------------------------------------------------
# 2. Magnetic feedback torque
# ---------------------------------------------------------------------------

def current_density(B_grid: np.ndarray,
                    spacing: float) -> np.ndarray:
    """
    Current density from Ampere's law, J = curl(B) / mu_0.

    Args:
        B_grid: magnetic field on a uniform 3D grid, shape (3, nx, ny, nz), T
        spacing: uniform grid spacing, m

    Returns:
        Current density on the same grid, shape (3, nx, ny, nz), A/m^2

    Raises:
        ValueError: if B_grid does not have leading dimension 3.
    """
    B_grid = np.asarray(B_grid, dtype=float)
    if B_grid.ndim != 4 or B_grid.shape[0] != 3:
        raise ValueError(
            f"B_grid must have shape (3, nx, ny, nz), got {B_grid.shape}")

    Bx, By, Bz = B_grid
    # np.gradient returns derivatives along each axis in order (x, y, z)
    dBx_dy = np.gradient(Bx, spacing, axis=1)
    dBx_dz = np.gradient(Bx, spacing, axis=2)
    dBy_dx = np.gradient(By, spacing, axis=0)
    dBy_dz = np.gradient(By, spacing, axis=2)
    dBz_dx = np.gradient(Bz, spacing, axis=0)
    dBz_dy = np.gradient(Bz, spacing, axis=1)

    curl = np.stack([
        dBz_dy - dBy_dz,
        dBx_dz - dBz_dx,
        dBy_dx - dBx_dy,
    ])
    return curl / MU_0


def lorentz_force_density(J: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Lorentz force density F = J x B.

    Args:
        J: current density, shape (3, ...) A/m^2
        B: magnetic field, shape (3, ...) T

    Returns:
        Force density, shape (3, ...), N/m^3
    """
    J = np.asarray(J, dtype=float)
    B = np.asarray(B, dtype=float)
    return np.cross(J, B, axis=0)


def magnetic_torque(r: np.ndarray,
                    J: np.ndarray,
                    B: np.ndarray,
                    cell_volume: float) -> np.ndarray:
    """
    Volume-integrated magnetic torque, tau = integral( r x (J x B) ) dV.

    Args:
        r: position vectors, shape (3, ...) m
        J: current density, shape (3, ...) A/m^2
        B: magnetic field, shape (3, ...) T
        cell_volume: volume of one grid cell, m^3

    Returns:
        Net torque vector, shape (3,), N*m
    """
    F = lorentz_force_density(J, B)
    torque_density = np.cross(np.asarray(r, dtype=float), F, axis=0)
    return torque_density.reshape(3, -1).sum(axis=1) * cell_volume


def effective_torque(lat_deg: float,
                     lon_deg: float,
                     B_magnitude: float,
                     v_magnitude: float,
                     beta: float = 1.0) -> float:
    """
    Roadmap's simplified hemispheric coupling torque:

        tau_eff = beta * sin(lat) * cos(lon) * |B| * |v|

    Antisymmetric in latitude, so it transfers coupling between hemispheres
    rather than adding net torque over a symmetric shell.

    Args:
        lat_deg: latitude, degrees (-90 to 90)
        lon_deg: longitude, degrees
        B_magnitude: |B|, T
        v_magnitude: |v|, m/s
        beta: coupling coefficient

    Returns:
        Effective torque in the roadmap's own units.
    """
    return (beta
            * np.sin(np.deg2rad(lat_deg))
            * np.cos(np.deg2rad(lon_deg))
            * B_magnitude
            * v_magnitude)


# ---------------------------------------------------------------------------
# 3. Earth-ionosphere dynamo
# ---------------------------------------------------------------------------

def ionospheric_dynamo_spec(v_atmo: Sequence[float],
                            B: Sequence[float],
                            omega: float = EARTH_OMEGA,
                            layer_height: float = 1.10e5,
                            conductivity: float = 1.0e-3) -> float:
    """
    Roadmap formula, implemented verbatim:

        E_dynamo = omega_earth * |v_atmo x B| * h * sigma

    Args:
        v_atmo: neutral wind velocity vector, m/s
        B: magnetic field vector, T
        omega: Earth angular speed, rad/s
        layer_height: dynamo layer thickness, m
        conductivity: height-averaged conductivity, S/m

    Returns:
        Value in the roadmap's own (non-SI-closed) units. The product
        |v x B| * sigma * h is a sheet current density (A/m); multiplying by
        omega gives A/(m*s), not an energy. Use
        `ionospheric_dynamo_power_density` for W/m^2.
    """
    E = motional_electric_field(v_atmo, B)
    return omega * float(np.linalg.norm(E)) * layer_height * conductivity


def ionospheric_dynamo_power_density(v_atmo: Sequence[float],
                                     B: Sequence[float],
                                     layer_height: float = 1.10e5,
                                     conductivity: float = 1.0e-3) -> float:
    """
    Dimensionally closed companion: height-integrated Joule dissipation,

        P = sigma * |v_atmo x B|^2 * h

    Args:
        v_atmo: neutral wind velocity vector, m/s
        B: magnetic field vector, T
        layer_height: dynamo layer thickness, m
        conductivity: height-averaged Pedersen conductivity, S/m

    Returns:
        Dissipated power per unit area, W/m^2
    """
    E = motional_electric_field(v_atmo, B)
    return conductivity * float(np.dot(E, E)) * layer_height


def dynamo_sheet_current(v_atmo: Sequence[float],
                         B: Sequence[float],
                         layer_height: float = 1.10e5,
                         conductivity: float = 1.0e-3) -> float:
    """
    Height-integrated dynamo current, K = sigma * |v x B| * h.

    Returns:
        Sheet current density, A/m
    """
    E = motional_electric_field(v_atmo, B)
    return conductivity * float(np.linalg.norm(E)) * layer_height


# ---------------------------------------------------------------------------
# 4. Spatial zoning
# ---------------------------------------------------------------------------

class Zone(str, Enum):
    """
    The four roadmap buckets.

    Caveat: POLAR, MIDLAT and EQUATOR partition the sphere by latitude and
    are mutually exclusive; OCEANIC is a surface type that overlaps all
    three. Treat OCEANIC as a coupled reservoir, not a fourth latitude band.
    """
    POLAR = "polar"
    MIDLAT = "midlat"
    EQUATOR = "equator"
    OCEANIC = "oceanic"


ZONES: Tuple[Zone, ...] = (Zone.POLAR, Zone.MIDLAT, Zone.EQUATOR, Zone.OCEANIC)

# Roadmap coupling matrix. Symmetric with a zero diagonal, so the exchange
# term below conserves total energy exactly.
COUPLING_MATRIX = np.array([
    [0.0, 0.3, 0.2, 0.1],   # Polar
    [0.3, 0.0, 0.4, 0.2],   # Midlat
    [0.2, 0.4, 0.0, 0.3],   # Equator
    [0.1, 0.2, 0.3, 0.0],   # Oceanic
])

# Fraction of Earth's surface in each latitude band, from the solid angle of
# the band: polar |lat|>60, midlat 30-60, equatorial |lat|<30. These three sum
# to 1. The oceanic entry is the ~71% ocean fraction and overlaps the others.
ZONE_AREA_FRACTION: Dict[Zone, float] = {
    Zone.POLAR: 1.0 - np.sin(np.deg2rad(60.0)),      # ~0.134
    Zone.MIDLAT: np.sin(np.deg2rad(60.0)) - np.sin(np.deg2rad(30.0)),  # ~0.366
    Zone.EQUATOR: np.sin(np.deg2rad(30.0)),          # 0.500
    Zone.OCEANIC: 0.71,
}


@dataclass
class SpatialCEED:
    """
    Zone-resolved CEED energy model.

    Replaces the scalar Earth of `convergence_model.py` with a vector over
    four zones coupled by `COUPLING_MATRIX`:

        dE_i/dt = source_i - sink_i * E_i + sum_j C_ij * (E_j - E_i)

    The exchange term is written as a difference, so it redistributes energy
    without creating or destroying any. With sources and sinks off, total
    energy is conserved to integration accuracy.
    """

    energy: np.ndarray = field(
        default_factory=lambda: np.array([100.0, 100.0, 100.0, 100.0]))
    coupling: np.ndarray = field(default_factory=lambda: COUPLING_MATRIX.copy())
    sources: np.ndarray = field(default_factory=lambda: np.zeros(4))
    sink_rates: np.ndarray = field(default_factory=lambda: np.zeros(4))

    def __post_init__(self):
        self.energy = np.asarray(self.energy, dtype=float)
        self.coupling = np.asarray(self.coupling, dtype=float)
        self.sources = np.asarray(self.sources, dtype=float)
        self.sink_rates = np.asarray(self.sink_rates, dtype=float)

        n = len(ZONES)
        if self.energy.shape != (n,):
            raise ValueError(f"energy must have shape ({n},), got {self.energy.shape}")
        if self.coupling.shape != (n, n):
            raise ValueError(f"coupling must have shape ({n},{n}), got {self.coupling.shape}")
        if np.any(np.diag(self.coupling) != 0.0):
            raise ValueError("coupling matrix must have a zero diagonal")

    def exchange(self, energy: np.ndarray) -> np.ndarray:
        """
        Net inter-zone energy flow, sum_j C_ij * (E_j - E_i).

        Returns:
            Per-zone net flow, same shape as `energy`. Sums to ~0 for a
            symmetric coupling matrix.
        """
        energy = np.asarray(energy, dtype=float)
        # differences[i, j] = E_j - E_i
        differences = energy[np.newaxis, :] - energy[:, np.newaxis]
        return np.sum(self.coupling * differences, axis=1)

    def derivative(self, energy: np.ndarray, t: float = 0.0) -> np.ndarray:
        """Per-zone energy derivative."""
        energy = np.asarray(energy, dtype=float)
        return self.sources - self.sink_rates * energy + self.exchange(energy)

    def step(self, dt: float) -> np.ndarray:
        """Advance one timestep with classical RK4."""
        E = self.energy
        k1 = self.derivative(E)
        k2 = self.derivative(E + 0.5 * dt * k1)
        k3 = self.derivative(E + 0.5 * dt * k2)
        k4 = self.derivative(E + dt * k3)
        self.energy = E + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        return self.energy

    def run(self, years: float = 3.0, dt: float = 1.0 / 12.0
            ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Integrate forward.

        Returns:
            t: time array, years
            history: shape (len(t), 4) zone energies over time
        """
        steps = int(round(years / dt))
        t = np.arange(steps + 1) * dt
        history = np.empty((steps + 1, len(ZONES)))
        history[0] = self.energy

        for i in range(1, steps + 1):
            history[i] = self.step(dt)

        return t, history

    def as_dict(self) -> Dict[str, float]:
        """Current energies keyed by zone name."""
        return {zone.value: float(e) for zone, e in zip(ZONES, self.energy)}


if __name__ == "__main__":
    sw = SolarWindBaseline()
    iono = IonosphereBaseline()

    print("CEED MHD / Spatial Module")
    print("=" * 60)

    print("\nSolar wind baseline (1 AU):")
    print(f"  speed              {sw.speed:.3e} m/s")
    print(f"  number density     {sw.number_density:.3e} m^-3")
    print(f"  mass density       {sw.mass_density:.3e} kg/m^3")
    print(f"  IMF magnitude      {sw.imf_magnitude:.3e} T")

    v = [sw.speed, 0.0, 0.0]
    B = [0.0, 0.0, sw.imf_magnitude]

    print("\nMHD energy injection:")
    print(f"  |v x B|            {np.linalg.norm(motional_electric_field(v, B)):.4f} V/m")
    print(f"  spec formula       {mhd_energy_injection(v, B, sw.mass_density):.4e} (roadmap units)")
    print(f"  Poynting flux      {poynting_flux(v, B):.4e} W/m^2")
    print(f"  kinetic flux       {solar_wind_kinetic_flux(sw.speed, sw.mass_density):.4e} W/m^2")

    v_atmo = [iono.neutral_wind_speed, 0.0, 0.0]
    B_iono = [0.0, 0.0, iono.field_magnitude]

    print("\nEarth-ionosphere dynamo:")
    print(f"  |v x B|            {np.linalg.norm(motional_electric_field(v_atmo, B_iono)):.4e} V/m")
    print(f"  spec formula       {ionospheric_dynamo_spec(v_atmo, B_iono):.4e} (roadmap units)")
    print(f"  sheet current      {dynamo_sheet_current(v_atmo, B_iono):.4e} A/m")
    print(f"  Joule dissipation  {ionospheric_dynamo_power_density(v_atmo, B_iono):.4e} W/m^2")

    print("\nEffective hemispheric torque (beta=1):")
    for lat in (-60, -30, 0, 30, 60):
        tau = effective_torque(lat, 0.0, iono.field_magnitude, iono.neutral_wind_speed)
        print(f"  lat {lat:+4d}          {tau:+.4e}")

    print("\nZone area fractions:")
    for zone, frac in ZONE_AREA_FRACTION.items():
        print(f"  {zone.value:9s}        {frac:.3f}")
    bands = sum(ZONE_AREA_FRACTION[z] for z in (Zone.POLAR, Zone.MIDLAT, Zone.EQUATOR))
    print(f"  (three latitude bands sum to {bands:.3f}; oceanic overlaps them)")

    print("\nSpatial run: energy injected at the pole only, 3 years")
    model = SpatialCEED(
        energy=np.array([200.0, 100.0, 100.0, 100.0]),
        sources=np.array([5.0, 0.0, 0.0, 0.0]),
        sink_rates=np.array([0.02, 0.02, 0.02, 0.02]),
    )
    t, history = model.run(years=3.0)
    print(f"  start  {np.array2string(history[0], precision=1)}  total {history[0].sum():.1f}")
    print(f"  end    {np.array2string(history[-1], precision=1)}  total {history[-1].sum():.1f}")

    print("\nConservation check: no sources, no sinks")
    closed = SpatialCEED(energy=np.array([300.0, 50.0, 100.0, 20.0]))
    t, hist = closed.run(years=10.0)
    print(f"  start  {np.array2string(hist[0], precision=2)}  total {hist[0].sum():.6f}")
    print(f"  end    {np.array2string(hist[-1], precision=2)}  total {hist[-1].sum():.6f}")
    print(f"  drift  {abs(hist[-1].sum() - hist[0].sum()):.2e}")
