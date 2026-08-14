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
    obs_ref=118.0,       # 2015 F10.7 annual mean (SC24 declining)
    E_ref=180.0,         # model initial energy
    scale=1.0,           # 1 sfu ≈ 1 energy unit (direct proxy)
    source='NOAA SWPC / NRCan Penticton Observatory',
)

MAGNETIC_SCALE = SubsystemScale(
    name='magnetic',
    observable='Kp geomagnetic index',
    units='Kp',
    obs_ref=1.9,         # 2015 annual mean Kp
    E_ref=92.5,          # model initial energy
    scale=0.05,          # Kp is compressed; 1 Kp unit ≈ 20 energy units
    source='GFZ Potsdam (https://kp.gfz.de/en/data)',
)

ATMOSPHERIC_SCALE = SubsystemScale(
    name='atmospheric',
    observable='Global mean surface temperature anomaly',
    units='K',
    obs_ref=1.09,        # 2015 anomaly above 1850-1900 (GISTEMP+0.19)
    E_ref=118.0,         # model initial energy
    scale=0.05,          # 1 K ≈ 20 energy units
    source='NASA GISTEMP v4 + 0.19K offset (IPCC AR6 1850-1900 base)',
)

OCEANIC_SCALE = SubsystemScale(
    name='oceanic',
    observable='Ocean heat content anomaly (0-700m)',
    units='10^22 J',
    obs_ref=14.5,        # 2015 OHC anomaly (relative to 1955-2006 mean)
    E_ref=110.0,         # model initial energy
    scale=0.3,           # 1 × 10^22 J ≈ 3.3 energy units
    source='NOAA/NCEI Levitus et al., Cheng et al. (2024)',
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
# Source: NOAA SWPC / NRCan Penticton Observatory
# SC24 peaked ~2014, minimum ~2019, SC25 rising through 2024
OBSERVED_F107 = {
    2010: 80,  2011: 113, 2012: 120, 2013: 123, 2014: 146,
    2015: 118, 2016: 89,  2017: 77,  2018: 70,  2019: 69,
    2020: 70,  2021: 88,  2022: 113, 2023: 152, 2024: 180,
}

# Ap geomagnetic index annual mean [nT] (linear equivalent of Kp)
# Source: GFZ Potsdam (https://kp.gfz.de/en/data)
# Ap is preferred for averaging because Kp is quasi-logarithmic.
# Approximate Kp equivalents shown in comments.
OBSERVED_AP = {
    2010: 5,  2011: 7,  2012: 9,  2013: 8,  2014: 9,   # Kp ~1.3-1.9
    2015: 9,  2016: 8,  2017: 8,  2018: 6,  2019: 5,   # Kp ~1.3-1.9
    2020: 5,  2021: 7,  2022: 10, 2023: 13, 2024: 14,   # Kp ~1.3-2.3
}

# Kp annual mean [Kp units]
#
# NOTE: these are annual means of the 3-hourly Kp values as published, NOT
# values converted from OBSERVED_AP above. A previous comment here claimed
# "Kp ≈ 0.3 * Ap^0.55"; that expression does not reproduce this table (it is
# low by up to 1.07 Kp units) and has been removed rather than left to imply
# a derivation that was never performed.
#
# Averaging Kp directly is itself a compromise: Kp is quasi-logarithmic, so
# the mean of Kp is not the Kp of the mean ap. OBSERVED_AP is retained above
# for anyone wanting to redo this properly via the official ap<->Kp scale.
OBSERVED_KP = {
    2010: 1.3, 2011: 1.7, 2012: 1.9, 2013: 1.8, 2014: 1.9,
    2015: 1.9, 2016: 1.8, 2017: 1.8, 2018: 1.5, 2019: 1.3,
    2020: 1.3, 2021: 1.7, 2022: 2.0, 2023: 2.3, 2024: 2.3,
}

# Global mean surface temperature anomaly [K above 1850-1900]
# Source: NASA GISTEMP v4 (1951-1980 base) + 0.19 K offset per IPCC AR6
# Cross-checked with HadCRUT5, Berkeley Earth
OBSERVED_TEMP = {
    2010: 0.91, 2011: 0.80, 2012: 0.84, 2013: 0.87, 2014: 0.94,
    2015: 1.09, 2016: 1.20, 2017: 1.11, 2018: 1.04, 2019: 1.17,
    2020: 1.20, 2021: 1.04, 2022: 1.08, 2023: 1.36, 2024: 1.47,
}

# Ocean heat content 0-700m [10^22 J, relative to 1955-2006 mean]
# Source: NOAA/NCEI Levitus et al., Cheng et al. (2024, 2025)
# Trend: ~1.0-1.5 x10^22 J/yr increase
OBSERVED_OHC = {
    2010: 9.5,  2011: 10.5, 2012: 11.0, 2013: 12.0, 2014: 13.0,
    2015: 14.5, 2016: 14.5, 2017: 14.0, 2018: 15.0, 2019: 16.0,
    2020: 17.5, 2021: 18.5, 2022: 19.5, 2023: 21.0, 2024: 22.5,
}

# Global primary energy consumption [EJ/yr]
# Source: Energy Institute Statistical Review (substitution method)
# Note: 2020 dip is COVID-19 pandemic effect
OBSERVED_ENERGY = {
    2010: 524, 2011: 531, 2012: 541, 2013: 550, 2014: 556,
    2015: 560, 2016: 566, 2017: 576, 2018: 590, 2019: 595,
    2020: 566, 2021: 604, 2022: 608, 2023: 620, 2024: 632,
}

ALL_OBSERVED = {
    'solar': OBSERVED_F107,
    'magnetic': OBSERVED_KP,
    'atmospheric': OBSERVED_TEMP,
    'oceanic': OBSERVED_OHC,
}


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


def energy_to_observations(system: str, E: np.ndarray) -> np.ndarray:
    """Convert model energy trajectory to physical units.

    Args:
        system: Subsystem name.
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
