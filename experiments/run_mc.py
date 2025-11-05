“””
Monte Carlo Uncertainty Analysis
Sample over parameter ranges to quantify outcome uncertainty
“””

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class ParameterRange:
“”“Parameter with uncertainty range”””
name: str
mean: float
low: float
high: float
distribution: str = ‘uniform’  # ‘uniform’ or ‘normal’

def sample_parameter(param: ParameterRange) -> float:
“”“Sample a single value from parameter distribution”””
if param.distribution == ‘uniform’:
return np.random.uniform(param.low, param.high)
elif param.distribution == ‘normal’:
# Use mean and interpret range as ±2σ
sigma = (param.high - param.low) / 4.0
return np.random.normal(param.mean, sigma)
else:
return param.mean

def run_monte_carlo(n_samples: int = 200, horizon_years: float = 10) -> dict:
“””
Run Monte Carlo simulation over parameter uncertainty

```
Args:
    n_samples: Number of Monte Carlo samples
    horizon_years: Simulation duration
    
Returns:
    Dictionary with results arrays
"""
# Define uncertain parameters (IPCC AR6 ranges)
params = {
    'ECS': ParameterRange('ECS', mean=3.0, low=2.5, high=4.0, distribution='normal'),
    'aerosol_ERF': ParameterRange('Aerosol ERF', mean=-1.1, low=-1.7, high=-0.4, distribution='uniform'),
    'sink_strength': ParameterRange('Sink strength', mean=0.54, low=0.45, high=0.60, distribution='uniform'),
    'permafrost_rate': ParameterRange('Permafrost', mean=95, low=14, high=175, distribution='uniform'),
}

# Storage for results
final_temps = []
final_co2 = []
parameter_samples = {name: [] for name in params.keys()}

print(f"Running {n_samples} Monte Carlo samples...")

for i in range(n_samples):
    if (i + 1) % 50 == 0:
        print(f"  Sample {i+1}/{n_samples}")
    
    # Sample parameters
    sampled_params = {}
    for name, param_range in params.items():
        value = sample_parameter(param_range)
        sampled_params[name] = value
        parameter_samples[name].append(value)
    
    # Run simplified model with sampled parameters
    T_final, CO2_final = run_simple_model(
        horizon_years=horizon_years,
        **sampled_params
    )
    
    final_temps.append(T_final)
    final_co2.append(CO2_final)

return {
    'temperatures': np.array(final_temps),
    'co2': np.array(final_co2),
    'parameters': parameter_samples,
    'n_samples': n_samples,
    'horizon': horizon_years
}
```

def run_simple_model(horizon_years: float, ECS: float, aerosol_ERF: float,
sink_strength: float, permafrost_rate: float) -> Tuple[float, float]:
“””
Simplified model for fast Monte Carlo sampling

```
Args:
    horizon_years: Simulation duration
    ECS: Equilibrium climate sensitivity
    aerosol_ERF: Aerosol forcing
    sink_strength: Carbon sink effectiveness
    permafrost_rate: Permafrost feedback strength
    
Returns:
    (final_temperature, final_CO2)
"""
# Initial conditions
T = 1.1  # Current temp above pre-industrial
CO2 = 420  # Current CO2 ppm
CO2_preindustrial = 280

# Timestep
dt = 0.1
steps = int(horizon_years / dt)

for _ in range(steps):
    # CO2 forcing
    CO2_forcing = ECS * np.log(CO2 / CO2_preindustrial) / np.log(2)
    
    # Permafrost feedback
    if T > 0.5:
        permafrost_forcing = (permafrost_rate / 1000.0) * 0.5 * (T - 0.5)
    else:
        permafrost_forcing = 0.0
    
    # Total forcing
    total_forcing = CO2_forcing + aerosol_ERF + permafrost_forcing
    
    # Retention (collapses with temperature)
    retention = 1.05 * np.exp(-0.0008 * T**2)
    
    # Dissipation
    dissipation = 0.04 * T + 0.01 * (T**1.2)
    
    # Temperature update
    dT = 0.1 * (total_forcing * retention - dissipation) * dt
    T += dT
    
    # CO2 update
    emissions = 10.0  # GtC/year
    sink = sink_strength * np.exp(-0.08 * T)  # Weakens with warming
    net_emissions = emissions * (1 - sink)
    dCO2 = (net_emissions / 2.12) * dt
    CO2 += dCO2

return T, CO2
```

def analyze_results(results: dict):
“””
Analyze and visualize Monte Carlo results

```
Args:
    results: Dictionary from run_monte_carlo
"""
temps = results['temperatures']
co2 = results['co2']

print("\n" + "=" * 60)
print("Monte Carlo Results")
print("=" * 60)

# Temperature statistics
print(f"\nTemperature after {results['horizon']} years:")
print(f"  Mean: {np.mean(temps):.2f}°C")
print(f"  Median: {np.median(temps):.2f}°C")
print(f"  5th percentile: {np.percentile(temps, 5):.2f}°C")
print(f"  95th percentile: {np.percentile(temps, 95):.2f}°C")
print(f"  Range: {np.min(temps):.2f} - {np.max(temps):.2f}°C")

# Threshold exceedance
exceed_15 = np.sum(temps > 1.5) / len(temps) * 100
exceed_20 = np.sum(temps > 2.0) / len(temps) * 100
exceed_30 = np.sum(temps > 3.0) / len(temps) * 100

print(f"\nProbability of exceeding:")
print(f"  1.5°C: {exceed_15:.1f}%")
print(f"  2.0°C: {exceed_20:.1f}%")
print(f"  3.0°C: {exceed_30:.1f}%")

# CO2 statistics
print(f"\nCO₂ after {results['horizon']} years:")
print(f"  Mean: {np.mean(co2):.1f} ppm")
print(f"  Range: {np.min(co2):.1f} - {np.max(co2):.1f} ppm")
```

def plot_results(results: dict, save_path: str = None):
“””
Create visualization of Monte Carlo results

```
Args:
    results: Dictionary from run_monte_carlo
    save_path: Optional path to save figure
"""
temps = results['temperatures']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Histogram
ax1.hist(temps, bins=30, alpha=0.7, edgecolor='black')
ax1.axvline(np.median(temps), color='red', linestyle='--', 
            label=f'Median: {np.median(temps):.2f}°C')
ax1.axvline(1.5, color='orange', linestyle='--', alpha=0.5, label='1.5°C target')
ax1.axvline(2.0, color='red', linestyle='--', alpha=0.5, label='2.0°C limit')
ax1.set_xlabel('Temperature (°C above pre-industrial)', fontsize=11)
ax1.set_ylabel('Frequency', fontsize=11)
ax1.set_title(f'Temperature Distribution (n={results["n_samples"]})', fontsize=12)
ax1.legend()
ax1.grid(alpha=0.3)

# Cumulative distribution
sorted_temps = np.sort(temps)
cdf = np.arange(1, len(sorted_temps) + 1) / len(sorted_temps)
ax2.plot(sorted_temps, cdf * 100, linewidth=2)
ax2.axvline(1.5, color='orange', linestyle='--', alpha=0.5, label='1.5°C')
ax2.axvline(2.0, color='red', linestyle='--', alpha=0.5, label='2.0°C')
ax2.set_xlabel('Temperature (°C above pre-industrial)', fontsize=11)
ax2.set_ylabel('Cumulative Probability (%)', fontsize=11)
ax2.set_title('Cumulative Distribution', fontsize=12)
ax2.legend()
ax2.grid(alpha=0.3)

plt.tight_layout()

if save_path:
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
else:
    plt.show()
```

def main():
import argparse

```
parser = argparse.ArgumentParser(description='Run Monte Carlo uncertainty analysis')
parser.add_argument('--n', type=int, default=200, help='Number of samples')
parser.add_argument('--horizon', type=float, default=10, help='Simulation horizon (years)')
parser.add_argument('--plot', action='store_true', help='Show plots')
parser.add_argument('--save', type=str, default=None, help='Save plot to file')

args = parser.parse_args()

# Run Monte Carlo
results = run_monte_carlo(n_samples=args.n, horizon_years=args.horizon)

# Analyze
analyze_results(results)

# Plot
if args.plot or args.save:
    plot_results(results, save_path=args.save)
```

if **name** == “**main**”:
main()
