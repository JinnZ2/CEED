"""
Minimal Earth System Model (ESM)
IPCC AR6-calibrated two-variable energy balance model.

State variables:
    T   : global mean surface temperature anomaly [deg C above pre-industrial]
    CO2 : atmospheric CO2 concentration [ppm]

Governing equations:

    C dT/dt = F_total(T, CO2, t) - lambda_eff * T

    dCO2/dt = E_net(T) / alpha_CO2

where:
    C           = effective heat capacity of climate system [W yr / (m^2 K)]
    F_total     = sum of radiative forcings [W/m^2]
    lambda_eff  = effective climate feedback parameter [W / (m^2 K)]
    E_net       = net CO2 emissions after sink uptake [GtC/yr]
    alpha_CO2   = airborne fraction conversion (2.12 GtC per ppm)

Forcings included:
    - CO2 radiative forcing (logarithmic, IPCC AR6 formula)
    - Solar cycle (11-year cosine modulation)
    - Aerosol effective radiative forcing (scenario-dependent)
    - Permafrost carbon feedback (threshold-activated)
    - Cloud feedback (saturating positive)

References:
    - IPCC AR6 WG1 Chapter 7 (energy budget, ECS, feedback parameter)
    - Myhre et al. (1998): CO2 radiative forcing formula
    - Forster et al. (2021): ERF assessment in AR6
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Tuple
import argparse


@dataclass
class ClimateParameters:
    """IPCC AR6-informed parameters with stated ranges.

    Default values are best estimates; ranges are for uncertainty analysis.
    """

    # Equilibrium Climate Sensitivity [deg C per CO2 doubling]
    # AR6 WG1 Ch7: likely 2.5-4.0, best estimate 3.0
    ECS: float = 3.0

    # Climate feedback parameter [W/(m^2 K)]
    # lambda = F_2xCO2 / ECS, where F_2xCO2 ~ 3.7 W/m^2
    # At ECS=3.0: lambda ~ 1.23 W/(m^2 K)
    @property
    def lambda_eff(self) -> float:
        F_2xCO2 = 3.7  # W/m^2, forcing from CO2 doubling (Myhre et al. 1998)
        return F_2xCO2 / self.ECS

    # Effective heat capacity [W yr / (m^2 K)]
    # Ocean mixed layer ~100m => C ~ 4.2e8 J/(m^2 K) ~ 13.3 W yr/(m^2 K)
    # Using smaller value for upper-ocean response timescale (~10 yr)
    heat_capacity: float = 10.0

    # Aerosol effective radiative forcing [W/m^2], 1750-2019
    # AR6: best estimate -1.1, very likely -1.7 to -0.4
    aerosol_ERF: float = -1.1

    # CO2 radiative forcing coefficient [W/m^2]
    # F = F_2xCO2 * ln(CO2/CO2_0) / ln(2)
    F_2xCO2: float = 3.7

    # Carbon sinks
    # AR6: land + ocean absorb ~54% of anthropogenic CO2
    total_sink_fraction: float = 0.54
    sink_weakening_rate: float = 0.08  # exponential decay with T [1/K]

    # Permafrost feedback [GtCO2 per K above 0.5 K threshold]
    # AR6: wide range 14-175 GtCO2/K, mid-estimate ~95
    permafrost_sensitivity: float = 95.0
    permafrost_threshold: float = 0.5  # K above pre-industrial

    # Solar cycle
    solar_cycle_period: float = 11.0  # years
    solar_amplitude: float = 0.1  # W/m^2 peak-to-peak

    # Cloud feedback [W/(m^2 K)]
    # AR6: net cloud feedback +0.45 W/(m^2 K), range -0.1 to +0.97
    cloud_feedback_strength: float = 0.45
    cloud_saturation_temp: float = 4.0  # K where cloud feedback saturates

    # Emissions [GtC/yr]
    # Global Carbon Budget 2024: total anthropogenic (fossil + LULUCF) was
    # 11.1 GtC/yr in 2023, preliminary 11.3 GtC/yr in 2024.
    # Was 10.0 with a comment reading "current ~11 GtC/yr" — the comment was
    # right and the value was not.
    emissions_rate: float = 11.1

    # Unit conversion
    GtC_per_ppm: float = 2.12  # 1 ppm CO2 ~ 2.12 GtC in atmosphere
    CO2_preindustrial: float = 280.0  # ppm


class MinimalESM:
    """Two-variable Earth System Model: (T, CO2).

    Uses the standard energy balance formulation:
        C dT/dt = F_net - lambda * T
    """

    def __init__(self, params: ClimateParameters = None):
        self.params = params or ClimateParameters()
        self.T_initial = 1.1    # deg C above pre-industrial (~2024)
        self.CO2_initial = 420  # ppm (~2024)
        self.aerosol_scenario = "current"

    def co2_forcing(self, CO2: float) -> float:
        """Radiative forcing from CO2 [W/m^2].

        F_CO2 = F_2xCO2 * ln(CO2 / CO2_0) / ln(2)

        Reference: Myhre et al. (1998), adopted by IPCC AR6.
        """
        p = self.params
        if CO2 <= 0:
            return 0.0
        return p.F_2xCO2 * np.log(CO2 / p.CO2_preindustrial) / np.log(2)

    def solar_forcing(self, t: float) -> float:
        """Solar cycle forcing [W/m^2].

        Cosine modulation with 11-year period.
        """
        p = self.params
        return p.solar_amplitude * np.cos(2 * np.pi * t / p.solar_cycle_period)

    def aerosol_forcing(self, t: float) -> float:
        """Aerosol forcing [W/m^2] under selected policy scenario."""
        p = self.params

        if self.aerosol_scenario == "current":
            return p.aerosol_ERF

        elif self.aerosol_scenario == "regulated":
            # Linear reduction from current to -0.5 W/m^2 over 10 years
            target = -0.5
            rate = (p.aerosol_ERF - target) / 10.0
            return max(target, p.aerosol_ERF - rate * t)

        elif self.aerosol_scenario == "removed":
            # Rapid removal over 3 years
            rate = p.aerosol_ERF / 3.0
            return min(0.0, p.aerosol_ERF - rate * t)

        return 0.0

    def cloud_feedback_forcing(self, T: float) -> float:
        """Cloud feedback contribution [W/m^2].

        Saturating positive feedback:
            F_cloud = alpha * T * (1 - T / T_sat)  for T in [0, T_sat]

        Physical basis: low clouds decrease with warming (positive feedback),
        but effect saturates at high T as cloud regime shifts complete.
        """
        p = self.params
        if T <= 0:
            return 0.0
        f = p.cloud_feedback_strength * T * (1.0 - T / p.cloud_saturation_temp)
        return max(0.0, f)

    def permafrost_forcing(self, T: float) -> float:
        """Forcing from permafrost carbon release [W/m^2].

        Activated above threshold temperature.
        Conversion: permafrost_sensitivity GtCO2/K -> W/m^2 via
        approximate relationship: 1000 GtCO2 ~ 0.5 W/m^2 sustained.
        """
        p = self.params
        if T <= p.permafrost_threshold:
            return 0.0
        co2_released = p.permafrost_sensitivity * (T - p.permafrost_threshold)
        return (co2_released / 1000.0) * 0.5

    def carbon_sink_fraction(self, T: float) -> float:
        """Fraction of emissions absorbed by land + ocean sinks.

        Decreases exponentially with warming (AR6 WG1 Ch5).

        Returns value in [0, 1].
        """
        p = self.params
        return p.total_sink_fraction * np.exp(-p.sink_weakening_rate * T)

    def derivatives(self, state: np.ndarray, t: float) -> np.ndarray:
        """ODE right-hand side.

        State = [T, CO2]

        dT/dt  = (1/C) * [F_total - lambda * T]
        dCO2/dt = E_net / alpha_CO2
        """
        T, CO2 = state
        p = self.params

        # Total radiative forcing [W/m^2]
        F_total = (self.co2_forcing(CO2)
                   + self.solar_forcing(t)
                   + self.aerosol_forcing(t)
                   + self.cloud_feedback_forcing(T)
                   + self.permafrost_forcing(T))

        # Energy balance: C dT/dt = F_total - lambda * T
        dT_dt = (F_total - p.lambda_eff * T) / p.heat_capacity

        # CO2 budget
        sink = self.carbon_sink_fraction(T)
        net_emissions = p.emissions_rate * (1.0 - sink)  # GtC/yr
        dCO2_dt = net_emissions / p.GtC_per_ppm  # ppm/yr

        return np.array([dT_dt, dCO2_dt])

    def simulate(self, years: float = 10, dt: float = 0.1
                 ) -> Tuple[np.ndarray, np.ndarray]:
        """Integrate the model forward.

        Returns:
            t: time array [years]
            solution: state array [[T, CO2], ...], shape (N, 2)
        """
        t = np.arange(0, years, dt)
        state0 = np.array([self.T_initial, self.CO2_initial])
        solution = odeint(self.derivatives, state0, t)
        return t, solution

    def plot_results(self, t: np.ndarray, solution: np.ndarray,
                     save_path: str = None):
        """Plot temperature and CO2 trajectories."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        T = solution[:, 0]
        CO2 = solution[:, 1]

        ax1.plot(t, T, 'r-', linewidth=2, label='Temperature anomaly')
        ax1.axhline(y=1.5, color='orange', linestyle='--', alpha=0.5,
                     label='Paris 1.5 K target')
        ax1.axhline(y=2.0, color='red', linestyle='--', alpha=0.5,
                     label='Paris 2.0 K limit')
        ax1.set_ylabel('Temperature anomaly (K)', fontsize=12)
        ax1.set_title(
            f'Minimal ESM: {self.aerosol_scenario} aerosol scenario',
            fontsize=14)
        ax1.legend()
        ax1.grid(alpha=0.3)

        ax2.plot(t, CO2, 'b-', linewidth=2, label='CO2 concentration')
        ax2.axhline(y=450, color='orange', linestyle='--', alpha=0.5,
                     label='450 ppm threshold')
        ax2.set_xlabel('Years from present', fontsize=12)
        ax2.set_ylabel('CO2 (ppm)', fontsize=12)
        ax2.legend()
        ax2.grid(alpha=0.3)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        else:
            plt.show()


def main():
    parser = argparse.ArgumentParser(description='Run Minimal ESM')
    parser.add_argument('--horizon', type=float, default=10,
                        help='Simulation horizon (years)')
    parser.add_argument('--aerosol',
                        choices=['current', 'regulated', 'removed'],
                        default='current', help='Aerosol policy scenario')
    parser.add_argument('--plot', action='store_true', help='Show plot')
    parser.add_argument('--save', type=str, default=None,
                        help='Save plot to file')

    args = parser.parse_args()

    model = MinimalESM()
    model.aerosol_scenario = args.aerosol

    print(f"Running Minimal ESM for {args.horizon} years")
    print(f"Aerosol scenario: {args.aerosol}")
    print(f"ECS: {model.params.ECS} K,  lambda: {model.params.lambda_eff:.2f} W/(m^2 K)")
    print(f"Heat capacity: {model.params.heat_capacity} W yr/(m^2 K)")
    print("=" * 60)

    t, solution = model.simulate(years=args.horizon)

    T_final = solution[-1, 0]
    CO2_final = solution[-1, 1]

    print(f"\nResults after {args.horizon} years:")
    print(f"  Temperature: {T_final:.2f} K above pre-industrial "
          f"({T_final - model.T_initial:+.2f} K change)")
    print(f"  CO2: {CO2_final:.1f} ppm")

    if T_final > 2.0:
        print(f"  EXCEEDS Paris 2.0 K limit")
    elif T_final > 1.5:
        print(f"  EXCEEDS Paris 1.5 K target")
    else:
        print(f"  Within Paris targets")

    if args.plot or args.save:
        model.plot_results(t, solution, save_path=args.save)


if __name__ == "__main__":
    main()
