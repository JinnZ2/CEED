# CLAUDE.md

Guidelines for AI assistants working on the CEED codebase.

## Project Overview

CEED (Cascading Energetic Event Disruption) is a Python simulation framework
for modelling cross-domain energy convergence across Earth-space systems.

**What CEED tests that other models don't**: existing Earth System Models
(CESM, GFDL, UKESM) treat solar, magnetic, atmospheric, and oceanic systems
largely in isolation.  CEED's thesis is that cross-domain energy transfer
pathways can amplify anthropogenic forcing beyond single-system predictions.
The model's value is in the **coupling matrix**, not in beating full ESMs at
single-system accuracy.

Parameters are calibrated against 2010-2024 observations via the hindcast
framework.  See `Docs/calibration-guide.md` for tuning workflow.

## Repository Structure

```
CEED/
├── CEED_universal_model.py          # Domain-agnostic feedback framework
├── dashboard_starter.py             # Streamlit UI prototype
├── conftest.py                      # Puts repo root on sys.path for tests
├── simulation/
│   ├── convergence_model.py         # Multi-subsystem coupled ODE model
│   ├── minimum_esm_code.py          # 2-variable Earth System Model (T, CO2)
│   ├── mhd_spatial_model.py         # MHD injection, torque, dynamo, zone coupling
│   ├── unit_bridge.py               # Energy indices <-> physical observables
│   ├── hindcast.py                  # Validation against 2010-2024 observations
│   ├── validation.py                # Held-out train/test split
│   ├── tipping.py                   # Bistability, hysteresis, commitment
│   └── cascade.py                   # Coupled tipping elements (cascades)
├── Tests/                           # Pytest suite (bare `pytest` collects all)
├── experiments/
│   └── run_mc.py                    # Monte Carlo uncertainty analysis
├── Data/
│   └── inputs.py                    # Mock data input layer
├── Docs/
│   ├── CEED-model-specs.md          # Full model specification
│   ├── calibration-guide.md         # How to tune parameters against data
│   └── numerical-audit.md           # Parameter/threshold/units audit, open findings
├── legacy/                          # Superseded versions + falsification record
├── requirements.txt                 # Python dependencies
├── references.md                    # Scientific literature citations
└── red-team-report.md               # Risk assessment report
```

## Quick Reference

| What | Command / Detail |
|------|-----------------|
| Language | Python 3.8+ |
| Install deps | `pip install -r requirements.txt` |
| Run tests | `pytest` |
| Run simulation | `python simulation/convergence_model.py` |
| Run ESM | `python simulation/minimum_esm_code.py --plot` |
| Run Monte Carlo | `python experiments/run_mc.py --n 200 --horizon 10` |
| Run hindcast | `python -m simulation.hindcast --compare` |
| Held-out validation | `python -m simulation.validation` |
| Tipping elements | `python simulation/tipping.py` |
| Tipping cascades | `python simulation/cascade.py` |
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

### Unit Bridge (`simulation/unit_bridge.py`)

Maps between model energy indices and physical observables:

| Subsystem   | Observable                | Units    | Scale              |
|-------------|---------------------------|----------|--------------------|
| Solar       | F10.7 radio flux          | sfu      | 1 sfu ≈ 1 E-unit   |
| Magnetic    | Kp geomagnetic index      | Kp       | 1 Kp ≈ 20 E-units  |
| Atmospheric | Global mean temp anomaly  | K        | 1 K ≈ 20 E-units   |
| Oceanic     | Ocean heat content 0-700m | 10²² J   | 1×10²² J ≈ 3.3 E   |

Reference year: 2015 (model t=0). Historical data from NOAA, NASA GISS,
GFZ Potsdam, HadCRUT5 for 2010-2024.

### Hindcast Framework (`simulation/hindcast.py`)

Validates model against historical observations using:
- RMSE, MAE, R², NRMSE, Bias
- Letter grades (A-F) based on R² and NRMSE
- Baseline vs anthropogenic comparison

Run: `python -m simulation.hindcast --compare`

Current hindcast scores (anthropogenic model, 2010-2024):

| System      | R²    | Grade | Status                                  |
|-------------|-------|-------|-----------------------------------------|
| Oceanic     | +0.63 | B     | In-sample (tuned on this window)        |
| Atmospheric | +0.47 | C     | In-sample (tuned on this window)        |
| Solar       | +0.36 | C     | Fixed: phase lag derived, not fitted    |
| Magnetic    | -0.21 | F     | Open — cannot track the cycle           |
| **Mean**    | **+0.31** |   | Read the rows, not the average          |

**These are in-sample figures. Do not quote them as predictive skill.**

### Held-out validation (the number that matters)

`python -m simulation.validation` splits the record: train 2010-2019, test
2020-2024, one continuous trajectory, no re-initialisation at 2020.

| System      | train R² | test R² | gap    |
|-------------|----------|---------|--------|
| solar       | -0.020   | **+0.555** | -0.575 |
| magnetic    | -1.193   | +0.057  | -1.249 |
| atmospheric | +0.316   | +0.042  | +0.274 |
| oceanic     | **+0.773** | **-2.705** | **+3.478** |
| **Mean**    | -0.031   | **-0.513** | +0.482 |

Out of sample the shipped parameters are **beaten by a flat line** (mean test
R² -0.513). Oceanic — the best-looking subsystem in-sample — is the worst
offender, +0.773 train to -2.705 test.

Solar is the one that generalises (test +0.555), and it is the one whose
constants were **derived** from the model's linearisation rather than fitted.
That is the lesson: derived constants travelled, fitted ones did not.

Caveat: the test window is 5 years against an 11-year solar cycle, so this
tests drift and extrapolation, not cycle physics.

Two cautions when quoting any of these:

1. **Mean R² is a weak summary** across four incommensurable subsystems.
2. **Never quote the full-record scores as skill.** Use the held-out table.

Magnetic is structurally unable to oscillate: its source is constant and the
solar coupling carries only ~2% of its input. See finding 18 in
`Docs/numerical-audit.md`. It is deliberately untuned — fitting it against the
window used to score it would produce a better number and no new knowledge.

See `Docs/calibration-guide.md` for tuning priorities and workflow.

### Earth System Model (`simulation/minimum_esm_code.py`)

Standard energy balance form (IPCC AR6 WG1 Ch7):

```
C dT/dt = F_total(T, CO2, t) - lambda_eff * T
dCO2/dt = E_net(T) / alpha_CO2
```

- C = 10.0 W yr/(m^2 K) — effective heat capacity
- lambda_eff = F_2xCO2 / ECS — climate feedback parameter
- F_2xCO2 = 3.7 W/m^2 — CO2 doubling forcing (Myhre et al. 1998)

### Tipping Elements (`simulation/tipping.py`)

A fast-slow bistable primitive. The convergence model has **one attractor** —
four initial conditions spanning 250 to 1501 total energy all land within 0.05
of each other after 200 years — so it cannot represent a system pushed into a
new state that persists. This module supplies what is missing:

```
dx/dt = (x - x^3 + F(t) + noise) / tau_fast     internal state
dh/dt = (x - h) / tau_slow                       observable response
```

Fold at |F| = 2/(3 sqrt 3) = 0.3849, x = +/-1/sqrt(3).

- **Bistability**: two stable branches while |F| < 0.3849
- **Hysteresis**: up-ramp switches at +0.399, down-ramp at -0.399
- **Critical slowing down**: recovery time 0.50 -> 12.50 approaching the fold
- **Commitment lag**: x crosses while h has moved 3.2% of its eventual change

The commitment lag is the point. `tau_slow >> tau_fast` means the system is
decided long before it looks decided.

The shelf/sheet pair maps onto the two variables directly: `x` is ice SHELF
integrity (floating, disintegrates in weeks to months), `h` is grounded ice
SHEET mass (responds over centuries). They couple because the shelf buttresses
the sheet. Shelf loss raises sea level by almost nothing directly — it is
already floating — and matters because of what it stops holding back. The
visible event is fast and nearly harmless; the consequence it commits to is
slow and large.

**Event-induced tipping** (`escape_probability`, `barrier_height`): a system
can tip while the mean forcing stays below the fold, because a single discrete
excursion clears the barrier. The barrier collapses from 0.250 at F=0 to
0.00003 at F=0.384, so "safely below threshold" buys less and less.

At F=0 — barrier at maximum, mean forcing nowhere near the fold — 27% of
400-year trials still tip. Motivated by atmospheric rivers: ~3% of the time,
but 40-80% of winter meltwater on peninsula shelves, with measured rain on
Thwaites of 30 mm in summer and 9 mm in winter, driving hydrofracture.

Tail shape matters, but only when the fold sits far from the typical event
size. At 3.85x it, a heavy tail tips ~2x more often than a thin tail of
identical mean AND variance; at 1.28x the ordering reverses (E11). Variance
alone does not tell you the risk.

Noise is pre-generated and passed in, never drawn inside the derivative.

**Not integrated** into `ConvergencePredictor`. Wiring it in changes the
model's character and is a maintainer decision.

### Tipping Cascades (`simulation/cascade.py`)

CEED is named for cascades and until now modelled none. N coupled elements:

```
dx_i/dt = (x_i - x_i^3 + F_i(t) + sum_j C_ij*phi_j + noise) / tau_fast_i
dh_i/dt = (x_i - h_i) / tau_slow_i
phi_j   = (h_j + 1)/2                    element j's tipped fraction
```

`C_ij` is j's effect on i. Positive destabilises, negative protects. The
matrix is not symmetric — Thwaites->Ross is not Ross->Thwaites.

**Coupling runs through the SLOW variable.** Influence only arrives as the
source element actually responds, which makes cascades delay-dependent. That
is the whole point: coupling through x would erase the effect that matters.
Coupling uses `phi` rather than raw `h` so an intact network exerts no
influence at t=0.

Illustrative Thwaites/Ross configuration shows:

- uncoupled, tipping Thwaites leaves Ross untouched
- coupled at 0.55, Ross tips ~65 units later with no forcing of its own
- domino threshold: 0.3853 minimum coupling for the cascade
- Ross carries ~95% of the consequence and arrives last
- a stabilising Ross->Thwaites link cannot help, because Ross needs
  tau_slow=400 to respond and Thwaites tips at t=4. Sign is not enough;
  a protective coupling slower than the collapse is worthless.

Not calibrated. Do not read predictions out of it.

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

Tests use **pytest** and live in `Tests/`. A bare `pytest` from the repository root collects the whole suite.

Run: `pytest -v`

Key test properties verified (151 tests):
- Model initialization and output shape
- Phase classification range (0-4), thresholds relative to baseline
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
7. **Run tests**: Always run `pytest` after
   modifying simulation code.
8. **Hindcast first**: After any parameter change, run
   `python -m simulation.hindcast --compare` and check that scores
   don't regress.  The hindcast is the ground truth.
9. **No over-engineering**: This is a research/simulation project. Keep
   abstractions minimal and code readable to scientists.
