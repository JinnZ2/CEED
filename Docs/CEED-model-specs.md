# CEED Model Specification

**Cascading Energetic Event Disruption**
**Multi-System Energy Convergence Framework**
Version: `v1.0.0`
Author: CEED Co-Creator
Date: July 2025

---

## Purpose

CEED models energy accumulation and cross-domain coupling across Earth-space
systems.  It evaluates whether positive feedbacks can overcome dissipation,
driving the coupled system through a sequence of phase transitions toward
nonlinear amplification or cascade.

---

## 1. Convergence Model (Multi-Subsystem ODE)

File: `simulation/convergence_model.py`

### 1.1 State Variables

Four coupled subsystems, each characterised by an energy-like index E_i(t):

| Subsystem   | Proxy observable          | Baseline value |
|-------------|---------------------------|----------------|
| Solar       | F10.7 index (sfu)         | 180.0          |
| Magnetic    | Kp-derived index          | 92.5           |
| Atmospheric | Thermospheric density idx | 118.0          |
| Oceanic     | Ocean heat content idx    | 110.0          |

### 1.2 Governing Equation

Each subsystem obeys an energy balance ODE:

```
dE_i/dt = S_i(t) + (alpha_i - lambda_i) * E_i - gamma_i * E_i^2 + sum_j(c_ij * E_j)
```

where:

| Symbol       | Meaning                               | Units           |
|--------------|---------------------------------------|-----------------|
| S_i(t)       | Source / input rate                    | energy / yr     |
| alpha_i      | Retention rate (energy self-reinforcement) | 1 / yr      |
| lambda_i     | Linear dissipation rate               | 1 / yr          |
| gamma_i      | Nonlinear dissipation coefficient     | 1 / (energy·yr) |
| c_ij         | Cross-system coupling coefficient     | 1 / yr          |

**Energy conservation (1st law)**: energy entering the system either stays
(retention, alpha) or leaves (dissipation, lambda + gamma*E).  The net rate
`(alpha - lambda)` determines whether a subsystem accumulates or loses energy.

**2nd law constraint**: `gamma_i > 0` ensures dissipation dominates at high E
(entropy production scales superlinearly with energy flux), preventing
unphysical unbounded growth.  Solutions always converge to a finite
equilibrium.

### 1.3 Source Functions

```
S_solar(t)       = 5.0 * (1 + 0.3 * cos(2*pi*t / 11))   [11-yr solar cycle]
S_magnetic(t)    = -2.0 * (1 + 0.1*t)                    [secular field decay]
S_atmospheric(t) = 3.0 * (1 + 0.05*t)                    [anthropogenic trend]
S_oceanic(t)     = 1.0 * (1 + 0.02*t)                    [anthropogenic trend]
```

### 1.4 Retention and Dissipation Parameters

| Subsystem   | alpha_i | lambda_i | alpha-lambda | gamma_i |
|-------------|---------|----------|--------------|---------|
| Solar       | 0.06    | 0.05     | +0.01        | 0.0002  |
| Magnetic    | 0.025   | 0.02     | +0.005       | 0.0002  |
| Atmospheric | 0.09    | 0.08     | +0.01        | 0.0002  |
| Oceanic     | 0.015   | 0.01     | +0.005       | 0.0002  |

All subsystems have `alpha > lambda` (net positive retention), meaning they
accumulate energy at low-to-moderate levels.  This reflects the current
physical state: Earth is gaining ~1 W/m² net energy.

Physical basis for retention (alpha):
- **Solar**: coronal magnetic confinement of plasma
- **Magnetic**: ring current self-sustaining via gradient-curvature drift
- **Atmospheric**: greenhouse trapping of outgoing longwave radiation
- **Oceanic**: thermal inertia of deep ocean heat storage

The quadratic term `gamma_i * E_i^2` represents enhanced dissipation at
high energy (Stefan-Boltzmann T^4 scaling, enhanced convection, increased
particle precipitation).  It guarantees bounded solutions.

### 1.5 Cross-System Coupling

Only physically motivated couplings are nonzero:

| From → To              | c_ij  | Physical mechanism                |
|-------------------------|-------|-----------------------------------|
| Solar → Magnetic        | 0.005 | Solar wind drives geomagnetic storms |
| Solar → Atmospheric     | 0.003 | EUV heating of thermosphere       |
| Magnetic → Atmospheric  | 0.002 | Joule heating from auroral currents |
| Atmospheric → Oceanic   | 0.001 | Air-sea heat flux                 |

All other c_ij = 0.

### 1.6 External Forcing (Extended Model)

The `ExtendedConvergencePredictor` adds:

- **Meteor/bolide events**: Poisson process, rate 0.5/yr, energy 1–15 units.
- **Satellite launches**: ~150/yr, 2.0 energy units each.
- **Volcanic dissipation**: Poisson, ~0.6/yr, removes 5–25 units.

Events are pre-generated before integration (deterministic RHS required by
ODE solver) and applied as Gaussian pulses with width ~ 1 day.

An additional saturating sink represents uncharacterised geophysical losses:

```
D_unk(E_tot) = mu * E_tot * max(0.1, 1 - E_tot / E_sat)
```

with mu = 0.15 /yr and E_sat = 800 energy units.

---

## 2. Phase Classification

Phases are defined by total system energy E_total = sum(E_i):

| Phase | Threshold | Regime                   |
|-------|-----------|--------------------------|
| 1     | < 120     | Baseline                 |
| 2     | >= 120    | System stress            |
| 3     | >= 150    | Cross-system coupling    |
| 4     | >= 200    | Nonlinear amplification  |
| 5     | >= 300    | Cascade / collapse       |

Phase transitions occur when coupling terms `c_ij * E_j` become significant
relative to the dissipation terms, allowing energy injected in one subsystem
to amplify in others.

---

## 3. Minimal Earth System Model (ESM)

File: `simulation/minimum_esm_code.py`

### 3.1 State Variables

| Variable | Meaning                                   | Initial | Units |
|----------|-------------------------------------------|---------|-------|
| T        | Global mean surface temperature anomaly   | 1.1     | K     |
| CO2      | Atmospheric CO2 concentration             | 420     | ppm   |

### 3.2 Governing Equations

Standard energy balance form (IPCC AR6 WG1 Ch7):

```
C dT/dt = F_total(T, CO2, t) - lambda_eff * T

dCO2/dt = E_net(T) / alpha_CO2
```

| Symbol      | Value  | Units          | Source            |
|-------------|--------|----------------|-------------------|
| C           | 10.0   | W yr/(m^2 K)   | Ocean mixed layer |
| lambda_eff  | F_2x/ECS | W/(m^2 K)   | AR6 WG1 Ch7      |
| F_2xCO2    | 3.7    | W/m^2          | Myhre et al. 1998 |
| ECS         | 3.0    | K              | AR6 best estimate |
| alpha_CO2   | 2.12   | GtC/ppm        | Standard          |

### 3.3 Forcing Components

**CO2 radiative forcing** (logarithmic, Myhre et al. 1998):
```
F_CO2 = F_2xCO2 * ln(CO2 / 280) / ln(2)
```

**Solar cycle** (11-year cosine):
```
F_solar(t) = 0.1 * cos(2*pi*t / 11)   [W/m^2]
```

**Aerosol ERF**: -1.1 W/m^2 baseline (AR6: very likely -1.7 to -0.4).
Three scenarios: current, regulated (ramp to -0.5 over 10 yr), removed (ramp
to 0 over 3 yr).

**Cloud feedback** (saturating positive):
```
F_cloud = 0.45 * T * (1 - T/4.0)   [W/m^2]   for T > 0
```
AR6: net cloud feedback +0.45 W/(m^2 K), range -0.1 to +0.97.

**Permafrost** (threshold-activated):
```
F_pf = 0.5 * permafrost_sensitivity * (T - 0.5) / 1000   [W/m^2]   for T > 0.5 K
```
AR6: 14–175 GtCO2/K, mid-estimate 95.

### 3.4 Carbon Sink

Fraction of emissions absorbed by land + ocean:
```
f_sink(T) = 0.54 * exp(-0.08 * T)
```

AR6 baseline: 54% uptake, weakening exponentially with warming.

Net CO2 rate:
```
dCO2/dt = emissions * (1 - f_sink) / 2.12   [ppm/yr]
```

---

## 4. Universal Framework

File: `CEED_universal_model.py`

Domain-agnostic model for any feedback-driven system.

### 4.1 ODE

```
dE/dt = F_ext(t) + sum_k f_k(E) - D(E)
```

**Feedback rate** for loop k with strength s_k:
```
f_k(E) = (+/-) s_k * E * sigma(E, E_sat_k)
```
where `sigma(E, E_sat) = 1 / (1 + (E/E_sat)^2)` is a dimensionless
saturation function.

**Dissipation**:
```
D(E) = alpha * E + beta * |E|^p
```
Default: alpha = 0.05, beta = 0.001, p = 1.5.

### 4.2 Stability Metric

```
S = (sum_k f_k(E) + D(E)) / D(E)
```

- S > 1: net amplification (unstable)
- S = 1: balanced
- S < 1: net dissipation (stable)

### 4.3 Buffer Capacity

Buffer B in [0, 1] depletes above warning threshold E_warn:
```
dB/dt = -b * (E / E_warn) * B    when E > E_warn
```
with b = 0.01.

---

## 5. Monte Carlo Uncertainty Analysis

File: `experiments/run_mc.py`

Samples over IPCC AR6 parameter ranges using the same energy balance
formulation as the ESM:

| Parameter        | Distribution | Range           | Source     |
|------------------|-------------|-----------------|------------|
| ECS              | Normal      | 2.5 – 4.0 K    | AR6 likely |
| Aerosol ERF      | Uniform     | -1.7 – -0.4    | AR6 v.likely |
| Sink fraction    | Uniform     | 0.45 – 0.60    | AR6 WG1 Ch5 |
| Permafrost rate  | Uniform     | 14 – 175 GtCO2/K | AR6      |

Default: 200 samples, 10-year horizon.

---

## Data Sources

| Source                       | Data type                | Status   |
|------------------------------|--------------------------|----------|
| NOAA SWPC                    | Solar cycle, F10.7       | Real     |
| NASA OMNIWeb                 | Solar wind, IMF          | Real     |
| ESA Swarm                    | Geomagnetic field        | Real     |
| Emmert et al. (2020)         | Thermospheric density    | Cited    |
| Rahmstorf et al. (2023)      | AMOC weakening           | Cited    |
| IPCC AR6 WG1 (2021)          | ECS, ERF, feedbacks      | Cited    |
| Myhre et al. (1998)          | CO2 forcing formula      | Cited    |

Note: the `Data/Inputs` directory currently contains mock data.
Real-time API integration is planned.

---

## Limitations

1. Convergence model uses normalised energy indices, not physical units.
   Cross-system coupling coefficients are order-of-magnitude estimates.
2. ESM is a two-variable (T, CO2) reduced model; it omits ocean
   circulation dynamics, ice sheet dynamics, and regional variability.
3. Cloud and permafrost feedbacks use simplified functional forms.
4. Monte Carlo assumes parameter independence (no covariance structure).
