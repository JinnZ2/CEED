"""
Unit Bridge: mapping between convergence model energy indices
and physical observables.

The convergence model tracks normalised energy indices E_i(t) that proxy
for physical quantities measured by real instruments.  This module provides
the scaling functions and reference data needed to:

1. Convert physical observations → model energy indices
2. Convert model energy indices → physical units
3. Supply historical data for hindcast validation

Scaling approach:

    E_i = (obs_i - obs_ref_i) / scale_i + E_ref_i

where:
    obs_i     : physical observation in native units
    obs_ref_i : reference observation (e.g., long-term mean)
    scale_i   : conversion factor [physical_units / energy_unit]
    E_ref_i   : model reference energy (initial condition)

The inverse:
    obs_i = (E_i - E_ref_i) * scale_i + obs_ref_i

Physical observables and their model proxies:

    | Subsystem   | Observable                 | Units        | Source          |
    |-------------|----------------------------|--------------|-----------------|
    | Solar       | F10.7 radio flux           | sfu          | NOAA SWPC       |
    | Magnetic    | Kp geomagnetic index       | Kp units     | GFZ Potsdam     |
    | Atmospheric | Global mean temp anomaly   | K            | HadCRUT5/GISS   |
    | Oceanic     | Ocean heat content 0-700m  | 10^22 J      | NOAA/Levitus    |
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional


@dataclass
class SubsystemScale:
    """Scaling between a physical observable and model energy index.

    Attributes:
        name: Human-readable label.
        observable: What is measured (e.g., 'F10.7 radio flux').
        units: Physical units of the observable (e.g., 'sfu').
        obs_ref: Reference observation value.
        E_ref: Corresponding model energy index value.
        scale: Conversion factor [units / energy_index_unit].
        source: Data source citation.
    """
    name: str
    observable: str
    units: str
    obs_ref: float
    E_ref: float
    scale: float
    source: str

    def obs_to_energy(self, obs: float) -> float:
        """Convert physical observation to model energy index."""
        return (obs - self.obs_ref) / self.scale + self.E_ref

    def energy_to_obs(self, E: float) -> float:
        """Convert model energy index to physical observation."""
        return (E - self.E_ref) * self.scale + self.obs_ref

    def obs_array_to_energy(self, obs: np.ndarray) -> np.ndarray:
        """Vectorised observation → energy conversion."""
        return (obs - self.obs_ref) / self.scale + self.E_ref

    def energy_array_to_obs(self, E: np.ndarray) -> np.ndarray:
        """Vectorised energy → observation conversion."""
        return (E - self.E_ref) * self.scale + self.obs_ref


# ── Default scaling definitions ──────────────────────────────────

SOLAR_SCALE = SubsystemScale(
    name='solar',
    observable='F10.7 radio flux',
    units='sfu',
    obs_ref=150.0,       # approximate Solar Cycle 24-25 mean
    E_ref=180.0,         # model initial energy
    scale=1.0,           # 1 sfu ≈ 1 energy unit (direct proxy)
    source='NOAA SWPC / NRCan Ottawa',
)

MAGNETIC_SCALE = SubsystemScale(
    name='magnetic',
    observable='Kp geomagnetic index',
    units='Kp',
    obs_ref=2.0,         # quiet-time Kp
    E_ref=92.5,          # model initial energy
    scale=0.05,          # Kp is compressed; 1 Kp unit ≈ 20 energy units
    source='GFZ Potsdam / NOAA SWPC',
)

ATMOSPHERIC_SCALE = SubsystemScale(
    name='atmospheric',
    observable='Global mean surface temperature anomaly',
    units='K',
    obs_ref=1.1,         # ~2024 anomaly above pre-industrial
    E_ref=118.0,         # model initial energy
    scale=0.05,          # 1 K ≈ 20 energy units
    source='HadCRUT5 / NASA GISS / NOAA GlobalTemp',
)

OCEANIC_SCALE = SubsystemScale(
    name='oceanic',
    observable='Ocean heat content anomaly (0-700m)',
    units='10^22 J',
    obs_ref=15.0,        # ~2024 OHC anomaly (relative to 1955-2006 mean)
    E_ref=110.0,         # model initial energy
    scale=0.3,           # 1 × 10^22 J ≈ 3.3 energy units
    source='NOAA/NCEI Levitus et al.',
)

SCALES = {
    'solar': SOLAR_SCALE,
    'magnetic': MAGNETIC_SCALE,
    'atmospheric': ATMOSPHERIC_SCALE,
    'oceanic': OCEANIC_SCALE,
}


# ── Historical observation data ──────────────────────────────────
#
# Annual mean values for hindcast validation.
# Year 0 in the model = 2015 (midpoint of available data).
# t_model = year_calendar - 2015
#
# Sources and notes:
#   F10.7: NOAA SWPC / NRCan Ottawa (adjusted to 1 AU)
#   Kp: GFZ Potsdam Kp index (annual mean of 3-hourly values)
#   Temperature: NASA GISS LOTI (Land-Ocean Temperature Index)
#   OHC: NOAA/NCEI 0-700m (Levitus et al., relative to 1955-2006 mean)
#   Energy: IEA / BP Statistical Review (primary energy in EJ/yr)

REFERENCE_YEAR = 2015  # model t=0 corresponds to this calendar year

# F10.7 annual mean [sfu]
OBSERVED_F107 = {
    2010: 80,  2011: 113, 2012: 120, 2013: 123, 2014: 146,
    2015: 111, 2016: 89,  2017: 77,  2018: 70,  2019: 70,
    2020: 70,  2021: 90,  2022: 130, 2023: 160, 2024: 190,
}

# Kp annual mean
OBSERVED_KP = {
    2010: 1.3, 2011: 1.7, 2012: 1.6, 2013: 1.5, 2014: 1.8,
    2015: 1.8, 2016: 1.4, 2017: 1.2, 2018: 1.1, 2019: 1.0,
    2020: 1.0, 2021: 1.4, 2022: 1.7, 2023: 2.0, 2024: 2.3,
}

# Global mean surface temperature anomaly [K above 1850-1900]
# Source: NASA GISS (approximate annual means)
OBSERVED_TEMP = {
    2010: 0.72, 2011: 0.61, 2012: 0.64, 2013: 0.68, 2014: 0.75,
    2015: 0.87, 2016: 1.02, 2017: 0.92, 2018: 0.85, 2019: 0.98,
    2020: 1.02, 2021: 0.85, 2022: 0.89, 2023: 1.17, 2024: 1.29,
}

# Ocean heat content 0-700m [10^22 J, relative to 1955-2006 mean]
# Source: NOAA/NCEI Levitus (approximate annual)
OBSERVED_OHC = {
    2010: 8.0,  2011: 8.5,  2012: 9.5,  2013: 10.0, 2014: 11.0,
    2015: 11.5, 2016: 12.0, 2017: 12.0, 2018: 12.5, 2019: 13.5,
    2020: 14.5, 2021: 15.0, 2022: 15.5, 2023: 16.5, 2024: 17.5,
}

# Global primary energy consumption [EJ/yr]
# Source: IEA / BP Statistical Review
OBSERVED_ENERGY = {
    2010: 524, 2011: 541, 2012: 549, 2013: 558, 2014: 563,
    2015: 567, 2016: 572, 2017: 582, 2018: 596, 2019: 601,
    2020: 575, 2021: 604, 2022: 618, 2023: 620, 2024: 630,
}

ALL_OBSERVED = {
    'solar': OBSERVED_F107,
    'magnetic': OBSERVED_KP,
    'atmospheric': OBSERVED_TEMP,
    'oceanic': OBSERVED_OHC,
}


def years_to_model_time(years: Dict[int, float]) -> Dict[float, float]:
    """Convert calendar year keys to model time (t=0 at REFERENCE_YEAR)."""
    return {float(yr - REFERENCE_YEAR): val for yr, val in years.items()}


def observations_to_energy(system: str,
                           obs: Optional[Dict[int, float]] = None
                           ) -> Tuple[np.ndarray, np.ndarray]:
    """Convert observed data to model energy indices.

    Args:
        system: Subsystem name ('solar', 'magnetic', 'atmospheric', 'oceanic').
        obs: Dict of {calendar_year: observation}. If None, uses built-in data.

    Returns:
        t_model: model time array
        E_model: energy index array
    """
    if obs is None:
        obs = ALL_OBSERVED[system]

    scale = SCALES[system]
    t_model = np.array([float(yr - REFERENCE_YEAR) for yr in sorted(obs.keys())])
    obs_vals = np.array([obs[yr] for yr in sorted(obs.keys())])
    E_model = scale.obs_array_to_energy(obs_vals)

    return t_model, E_model


def energy_to_observations(system: str,
                           t_model: np.ndarray,
                           E: np.ndarray) -> np.ndarray:
    """Convert model energy trajectory to physical units.

    Args:
        system: Subsystem name.
        t_model: Model time array.
        E: Energy index array.

    Returns:
        Array of physical observations in native units.
    """
    return SCALES[system].energy_array_to_obs(E)


def print_scaling_table():
    """Print the unit bridge scaling table for documentation."""
    print("UNIT BRIDGE: Model Energy Index ↔ Physical Observable")
    print("=" * 75)
    print(f"{'System':<14s} {'Observable':<35s} {'Units':<10s} "
          f"{'Obs_ref':>8s} {'E_ref':>6s} {'Scale':>6s}")
    print("-" * 75)
    for name, s in SCALES.items():
        print(f"{name:<14s} {s.observable:<35s} {s.units:<10s} "
              f"{s.obs_ref:>8.1f} {s.E_ref:>6.1f} {s.scale:>6.2f}")
    print()
    print(f"Reference year: {REFERENCE_YEAR} (model t=0)")
    print(f"Conversion: E = (obs - obs_ref) / scale + E_ref")
    print(f"Inverse:    obs = (E - E_ref) * scale + obs_ref")


if __name__ == "__main__":
    print_scaling_table()

    print("\n\nSample conversions:")
    for name, s in SCALES.items():
        obs_2024 = list(ALL_OBSERVED[name].values())[-1]
        E = s.obs_to_energy(obs_2024)
        obs_back = s.energy_to_obs(E)
        print(f"  {name}: obs={obs_2024} {s.units} -> E={E:.1f} -> "
              f"obs={obs_back:.1f} {s.units}")
