# CLAUDE.md

Guidelines for AI assistants working on the CEED codebase.

## Project Overview

CEED (Cascading Energetic Event Disruption) is a Python simulation framework
for modelling cross-domain energy convergence across Earth-space systems.
It uses coupled ODEs to track energy accumulation in solar, magnetic,
atmospheric, and oceanic subsystems and evaluates whether positive feedbacks
can drive the system through phase transitions toward nonlinear amplification.

All parameters are anchored to IPCC AR6 and peer-reviewed literature.

## Repository Structure

```
CEED/
├── CEED_universal_model.py          # Domain-agnostic feedback framework
├── dashboard_starter.py             # Streamlit UI prototype
├── simulation/
│   ├── convergence_model.py         # Multi-subsystem coupled ODE model
│   └── minimum_esm_code.py          # 2-variable Earth System Model (T, CO2)
├── Tests/
│   └── test-convergence-model.py    # Pytest test suite
├── experiments/
│   └── run_mc.py                    # Monte Carlo uncertainty analysis
├── Data/
│   └── Inputs                       # Mock data input layer
├── Docs/
│   └── CEED-model-specs.md          # Full model specification
├── requirements.txt                 # Python dependencies
├── references.md                    # Scientific literature citations
└── red-team-report.md               # Risk assessment report
```

## Quick Reference

| What | Command / Detail |
|------|-----------------|
| Language | Python 3.8+ |
| Install deps | `pip install -r requirements.txt` |
| Run tests | `pytest Tests/test-convergence-model.py` |
| Run simulation | `python simulation/convergence_model.py` |
| Run ESM | `python simulation/minimum_esm_code.py --plot` |
| Run Monte Carlo | `python experiments/run_mc.py --n 200 --horizon 10` |
| Run dashboard | `streamlit run dashboard_starter.py` |
| License | MIT |

## Dependencies

Core: `numpy>=1.21.0`, `scipy>=1.7.0`, `matplotlib>=3.4.0`, `pyyaml>=5.4.0`

Optional: `streamlit` (dashboard), `pytest` (testing)

## Mathematical Formulations

### Convergence Model (`simulation/convergence_model.py`)

Four coupled subsystems (solar, magnetic, atmospheric, oceanic), each with
energy index E_i(t). The ODE is:

```
dE_i/dt = S_i(t) + A_i(t)
          + (alpha_i(E,t) - lambda_i) * E_i
          - gamma_i * E_i^2
          + sum_j [eta_ij(E) * c_ij(E) * E_j]   (received)
          - sum_j [c_ji(E) * E_i]                (sent)
```

- `S_i(t)`: natural source rate (solar cycle, secular trends)
- `A_i(t)`: anthropogenic forcing — geological energy release (fossil fuels)
- `alpha_i(E,t) * E_i`: retention — state-dependent, increases with
  anthropogenic load (more CO2 → more greenhouse trapping)
- `lambda_i * E_i`: linear dissipation (radiative loss)
- `gamma_i * E_i^2`: nonlinear dissipation — 2nd law bound
- `eta_ij * c_ij * E_j`: energy received (with conversion loss)
- `c_ji * E_i`: energy sent (full amount leaves)

**Anthropogenic forcing** (`AnthropogenicForcing`): models ~300 Myr of stored
solar energy released in ~200 yr. Logistic growth with resource peak.
Modifies not just source terms but system parameters themselves:

- Atmospheric retention increases logarithmically with load (greenhouse)
- Oceanic retention increases with surface warming (stratification)
- Coupling rates strengthen with energy gradients (Clausius-Clapeyron)
- Conversion efficiencies shift with total system energy

Coupling is **conservative**: 9 pathways with conversion efficiencies 5-40%.
Every transfer produces waste heat. See `Docs/CEED-model-specs.md` for
the full coupling matrix.

1st law: energy in = energy retained + energy dissipated + energy transferred.
When `alpha > lambda`, the system accumulates until `gamma * E^2` restores balance.

The extended model adds pre-generated external events (Gaussian pulses) and
a saturating unknown-sink term.

### Earth System Model (`simulation/minimum_esm_code.py`)

Standard energy balance form (IPCC AR6 WG1 Ch7):

```
C dT/dt = F_total(T, CO2, t) - lambda_eff * T
dCO2/dt = E_net(T) / alpha_CO2
```

- C = 10.0 W yr/(m^2 K) — effective heat capacity
- lambda_eff = F_2xCO2 / ECS — climate feedback parameter
- F_2xCO2 = 3.7 W/m^2 — CO2 doubling forcing (Myhre et al. 1998)

### Universal Framework (`CEED_universal_model.py`)

Domain-agnostic model:

```
dE/dt = F_ext(t) + sum_k f_k(E) - D(E)
```

Feedbacks use rational saturation: `f_k(E) = s_k * E / (1 + (E/E_sat)^2)`

Dissipation: `D(E) = alpha * E + beta * |E|^p`

## Code Style

No formatter or linter is configured. Follow these conventions:

- **Classes**: PascalCase (`ConvergencePredictor`, `SystemState`)
- **Functions/methods**: snake_case (`predict_convergence`, `classify_phases`)
- **Constants**: UPPER_CASE (`SYSTEMS`, `N_SYSTEMS`)
- **Type hints**: Use `typing` module and dataclass annotations
- **Data structures**: `@dataclass` for structured data
- **Imports**: stdlib, then third-party, grouped at file top
- Follow PEP 8

## Testing

Tests use **pytest** and live in `Tests/test-convergence-model.py`.

Run: `pytest Tests/test-convergence-model.py -v`

Key test properties verified (21 tests):
- Model initialization and output shape
- Phase classification range (1-4)
- ODE stays finite over long horizons (no blowup)
- Dissipation bounds growth (no source -> energy decays)
- ODE RHS is deterministic (no stochastic calls inside derivatives)
- Coupling is conservative (all eta in (0,1], waste heat produced)
- Extended model without events matches baseline exactly
- Anthropogenic forcing grows, peaks (logistic), increases total energy
- State-dependent params: no-anthro returns baseline, anthro increases
  retention and coupling, efficiency stays clamped [0.01, 0.95]

## Guidelines for AI Assistants

1. **ODE correctness**: The RHS of every ODE must be deterministic — no
   `random()` calls inside derivative functions. Pre-generate stochastic
   events before integration.
2. **Dimensional consistency**: Every term in `dE/dt` must have units of
   [energy/time]. Document units in docstrings and parameter tables.
3. **Thermodynamic consistency**: 1st law — retention (alpha) and
   dissipation (lambda) must both be present; energy doesn't vanish.
   2nd law — nonlinear dissipation (gamma*E^2) must dominate at high E,
   ensuring bounded solutions. Never remove retention without physical
   justification.
4. **State-dependent parameters**: When anthropogenic forcing is enabled,
   alpha, c, and eta become functions of (E, t). Any new parameter
   modification must have a physical mechanism documented in the docstring.
5. **Scientific accuracy**: Parameters must be traceable to IPCC AR6 or
   cited literature. Do not invent physical constants.
5. **Phase thresholds**: The 4-phase classification (120/150/200/300) is a
   core design choice. Do not change without explicit request.
6. **Standard forms**: Use `C dT/dt = F - lambda*T` for energy balance, not
   ad-hoc retention/dissipation multipliers.
7. **Run tests**: Always run `pytest Tests/test-convergence-model.py` after
   modifying simulation code.
8. **No over-engineering**: This is a research/simulation project. Keep
   abstractions minimal and code readable to scientists.
