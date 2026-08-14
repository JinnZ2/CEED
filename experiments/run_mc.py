"""
Monte Carlo Uncertainty Analysis for the Minimal ESM.

Samples over IPCC AR6 parameter ranges to quantify outcome uncertainty.
Uses the same energy balance formulation as minimum_esm_code.py:

    C dT/dt = F_total - lambda * T
    dCO2/dt = E_net / alpha_CO2

Parameters sampled:
    ECS             : Equilibrium Climate Sensitivity [K per CO2 doubling]
    aerosol_ERF     : Aerosol effective radiative forcing [W/m^2]
    sink_fraction   : Baseline carbon sink fraction [dimensionless]
    permafrost_rate : Permafrost carbon sensitivity [GtCO2/K]
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Tuple
import argparse


@dataclass
class ParameterRange:
    """Uncertain parameter with sampling distribution.

    Attributes:
        name: Human-readable label.
        mean: Central / best estimate.
        low: Lower bound of range.
        high: Upper bound of range.
        distribution: 'uniform' or 'normal'.
            For 'normal', range is interpreted as +/- 2 sigma.
    """
    name: str
    mean: float
    low: float
    high: float
    distribution: str = 'uniform'


def sample_parameter(param: ParameterRange) -> float:
    """Draw a single sample from the parameter distribution."""
    if param.distribution == 'uniform':
        return np.random.uniform(param.low, param.high)
    elif param.distribution == 'normal':
        sigma = (param.high - param.low) / 4.0
        return np.random.normal(param.mean, sigma)
    return param.mean


# Physical constants shared with minimum_esm_code.py
F_2xCO2 = 3.7          # W/m^2, forcing per CO2 doubling
CO2_PREINDUSTRIAL = 280 # ppm
GTC_PER_PPM = 2.12      # GtC per ppm CO2
HEAT_CAPACITY = 10.0    # W yr / (m^2 K)
EMISSIONS_RATE = 10.0   # GtC/yr


def run_single(horizon_years: float, ECS: float, aerosol_ERF: float,
               sink_fraction: float, permafrost_rate: float
               ) -> Tuple[float, float]:
    """Run a single forward integration with given parameters.

    Uses the same physics as MinimalESM but inlined for speed.

    Returns:
        (final_temperature, final_CO2)
    """
    T = 1.1    # current anomaly [K]
    CO2 = 420  # current CO2 [ppm]

    # Derived feedback parameter
    lambda_eff = F_2xCO2 / ECS  # W/(m^2 K)

    dt = 0.1   # years
    steps = int(horizon_years / dt)

    for _ in range(steps):
        # Forcings [W/m^2]
        F_co2 = F_2xCO2 * np.log(CO2 / CO2_PREINDUSTRIAL) / np.log(2)

        F_permafrost = 0.0
        if T > 0.5:
            F_permafrost = (permafrost_rate * (T - 0.5) / 1000.0) * 0.5

        # Cloud feedback (saturating, same as ESM)
        cloud_sat_T = 4.0
        F_cloud = 0.0
        if T > 0:
            F_cloud = max(0.0, 0.45 * T * (1.0 - T / cloud_sat_T))

        F_total = F_co2 + aerosol_ERF + F_permafrost + F_cloud

        # Energy balance: C dT/dt = F_total - lambda * T
        dT = (F_total - lambda_eff * T) / HEAT_CAPACITY * dt
        T += dT

        # CO2 budget
        sink = sink_fraction * np.exp(-0.08 * T)
        net_emissions = EMISSIONS_RATE * (1.0 - sink)
        dCO2 = (net_emissions / GTC_PER_PPM) * dt
        CO2 += dCO2

    return T, CO2


def run_monte_carlo(n_samples: int = 200,
                    horizon_years: float = 10) -> dict:
    """Run Monte Carlo ensemble over uncertain parameters.

    Returns dict with arrays of final temperatures and CO2.
    """
    params = {
        'ECS': ParameterRange(
            'ECS', mean=3.0, low=2.5, high=4.0, distribution='normal'),
        'aerosol_ERF': ParameterRange(
            'Aerosol ERF', mean=-1.1, low=-1.7, high=-0.4,
            distribution='uniform'),
        'sink_fraction': ParameterRange(
            'Sink fraction', mean=0.54, low=0.45, high=0.60,
            distribution='uniform'),
        'permafrost_rate': ParameterRange(
            'Permafrost', mean=95, low=14, high=175,
            distribution='uniform'),
    }

    final_temps = []
    final_co2 = []
    parameter_samples = {name: [] for name in params}

    print(f"Running {n_samples} Monte Carlo samples...")

    for i in range(n_samples):
        if (i + 1) % 50 == 0:
            print(f"  Sample {i + 1}/{n_samples}")

        sampled = {}
        for name, param_range in params.items():
            value = sample_parameter(param_range)
            sampled[name] = value
            parameter_samples[name].append(value)

        T_final, CO2_final = run_single(
            horizon_years=horizon_years,
            ECS=sampled['ECS'],
            aerosol_ERF=sampled['aerosol_ERF'],
            sink_fraction=sampled['sink_fraction'],
            permafrost_rate=sampled['permafrost_rate'],
        )

        final_temps.append(T_final)
        final_co2.append(CO2_final)

    return {
        'temperatures': np.array(final_temps),
        'co2': np.array(final_co2),
        'parameters': parameter_samples,
        'n_samples': n_samples,
        'horizon': horizon_years,
    }


def analyze_results(results: dict):
    """Print summary statistics for Monte Carlo ensemble."""
    temps = results['temperatures']

    print("\n" + "=" * 60)
    print("Monte Carlo Results")
    print("=" * 60)

    print(f"\nTemperature after {results['horizon']} years:")
    print(f"  Mean:             {np.mean(temps):.2f} K")
    print(f"  Median:           {np.median(temps):.2f} K")
    print(f"  5th percentile:   {np.percentile(temps, 5):.2f} K")
    print(f"  95th percentile:  {np.percentile(temps, 95):.2f} K")
    print(f"  Range:            {np.min(temps):.2f} - {np.max(temps):.2f} K")

    exceed_15 = np.mean(temps > 1.5) * 100
    exceed_20 = np.mean(temps > 2.0) * 100
    exceed_30 = np.mean(temps > 3.0) * 100

    print(f"\nExceedance probabilities:")
    print(f"  P(T > 1.5 K): {exceed_15:.1f}%")
    print(f"  P(T > 2.0 K): {exceed_20:.1f}%")
    print(f"  P(T > 3.0 K): {exceed_30:.1f}%")

    co2 = results['co2']
    print(f"\nCO2 after {results['horizon']} years:")
    print(f"  Mean:  {np.mean(co2):.1f} ppm")
    print(f"  Range: {np.min(co2):.1f} - {np.max(co2):.1f} ppm")


def plot_results(results: dict, save_path: str = None):
    """Histogram and CDF of final temperature distribution."""
    temps = results['temperatures']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.hist(temps, bins=30, alpha=0.7, edgecolor='black')
    ax1.axvline(np.median(temps), color='red', linestyle='--',
                label=f'Median: {np.median(temps):.2f} K')
    ax1.axvline(1.5, color='orange', linestyle='--', alpha=0.5,
                label='1.5 K target')
    ax1.axvline(2.0, color='red', linestyle='--', alpha=0.5,
                label='2.0 K limit')
    ax1.set_xlabel('Temperature (K above pre-industrial)', fontsize=11)
    ax1.set_ylabel('Frequency', fontsize=11)
    ax1.set_title(f'Temperature distribution (n={results["n_samples"]})',
                  fontsize=12)
    ax1.legend()
    ax1.grid(alpha=0.3)

    sorted_temps = np.sort(temps)
    cdf = np.arange(1, len(sorted_temps) + 1) / len(sorted_temps)
    ax2.plot(sorted_temps, cdf * 100, linewidth=2)
    ax2.axvline(1.5, color='orange', linestyle='--', alpha=0.5,
                label='1.5 K')
    ax2.axvline(2.0, color='red', linestyle='--', alpha=0.5,
                label='2.0 K')
    ax2.set_xlabel('Temperature (K above pre-industrial)', fontsize=11)
    ax2.set_ylabel('Cumulative probability (%)', fontsize=11)
    ax2.set_title('CDF', fontsize=12)
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Monte Carlo uncertainty analysis for Minimal ESM')
    parser.add_argument('--n', type=int, default=200,
                        help='Number of samples')
    parser.add_argument('--horizon', type=float, default=10,
                        help='Simulation horizon (years)')
    parser.add_argument('--plot', action='store_true', help='Show plots')
    parser.add_argument('--save', type=str, default=None,
                        help='Save plot to file')

    args = parser.parse_args()

    results = run_monte_carlo(n_samples=args.n, horizon_years=args.horizon)
    analyze_results(results)

    if args.plot or args.save:
        plot_results(results, save_path=args.save)


if __name__ == "__main__":
    main()
