# CLAUDE.md

Guidelines for AI assistants working on the CEED (Cascading Energetic Event Disruption) codebase.

## Project Overview

CEED is a Python simulation framework for modeling cross-domain energy convergence and feedback saturation across Earth-space systems. It explores worst-case scenarios where energy accumulation and positive feedback loops lead to nonlinear state transitions in climate, geomagnetic, oceanic, and atmospheric systems.

All parameters are anchored to IPCC AR6 Working Group I and peer-reviewed literature.

## Repository Structure

```
CEED/
├── CEED_universal_model.py          # Abstract feedback framework (FeedbackLoop, SystemState, CEEDSystem)
├── dashboard_starter.py             # Streamlit UI prototype
├── simulation/
│   ├── convergence_model.py         # Main extended model (ExtendedConvergencePredictor)
│   └── minimum_esm_code.py          # Minimal Earth System Model (IPCC-calibrated)
├── Tests/
│   └── test-convergence-model.py    # Pytest test suite
├── experiments/
│   └── run_mc.py                    # Monte Carlo uncertainty analysis
├── Data/
│   └── Inputs                       # Mock data input layer
├── Docs/
│   └── CEED-model-specs.md          # Model specification document
├── UI/                              # Placeholder for UI requirements
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
| Run dashboard | `streamlit run dashboard_starter.py` |
| Run Monte Carlo | `python experiments/run_mc.py` |
| License | MIT |

## Dependencies

Core: `numpy>=1.21.0`, `scipy>=1.7.0`, `matplotlib>=3.4.0`, `pyyaml>=5.4.0`

Optional: `streamlit` (dashboard UI, not in requirements.txt)

## Architecture

### Core Classes and Patterns

- **`FeedbackLoop`** (dataclass): Defines a feedback loop with `name`, `polarity`, `strength`, `saturation_threshold`.
- **`SystemState`** (dataclass): Tracks `energy`, `retention`, `dissipation`, `buffer_capacity`.
- **`CEEDSystem`**: Core engine that composes feedback loops and computes system evolution. Key metric: `stability_metric() = retention / dissipation` (>1.0 = unstable).
- **`ExtendedConvergencePredictor`**: Extended model tracking 4 subsystems (solar, magnetic, atmospheric, oceanic) with external event handling and phase classification.
- **Earth System Model** (`minimum_esm_code.py`): IPCC AR6-calibrated 2-layer energy balance model using `scipy.integrate.odeint`.

### Phase Classification

Energy thresholds define system phases (hardcoded):
- Phase 1: >= 120 (System Stress)
- Phase 2: >= 150 (Cross-System Coupling)
- Phase 3: >= 200 (Nonlinear Amplification)
- Phase 4: >= 300 (Cascade/Collapse)

### Energy Variables

- Energy variables follow the `E_*` naming pattern (e.g., `E_current`, `E_total`)
- Dissipation saturates at ~800 energy units
- Unknown sink baseline: 15% per timestep at low energy
- Retention collapse at high energy prevents unrealistic runaway

## Code Style

No formatter or linter is configured. Follow these observed conventions:

- **Classes**: PascalCase (`ConvergencePredictor`, `SystemState`)
- **Functions/methods**: snake_case (`predict_convergence`, `classify_phases`)
- **Constants**: UPPER_CASE
- **Type hints**: Use `typing` module (`List`, `Dict`, `Tuple`, `Callable`)
- **Data structures**: Prefer `@dataclass` for structured data
- **Docstrings**: Triple-quoted strings on classes and public methods
- **Imports**: Standard library first, then third-party, grouped at file top
- Follow PEP 8 conventions

## Testing

Tests use **pytest** and live in `Tests/test-convergence-model.py`. Current tests:

- `test_model_initialization` - Validates model initialization
- `test_energy_prediction_shape` - Checks prediction output dimensions
- `test_phase_classification` - Validates phase classification (1-4)
- `test_total_energy_increase` - Validates energy accumulation behavior

Run with: `pytest Tests/test-convergence-model.py`

There is no CI/CD pipeline configured. Always run tests locally before committing.

## Key Guidelines for AI Assistants

1. **Scientific accuracy**: All model parameters must be traceable to peer-reviewed literature or IPCC AR6. Do not invent physical constants.
2. **Preserve phase thresholds**: The 4-phase classification (120/150/200/300) is a core design choice. Do not change without explicit request.
3. **Energy conservation**: Models must respect energy balance. Dissipation saturation and retention collapse are intentional safeguards.
4. **Run tests**: Always run `pytest Tests/test-convergence-model.py` after modifying simulation code.
5. **No over-engineering**: This is a research/simulation project. Keep abstractions minimal and code readable to scientists.
6. **Data files**: Real data sources include NOAA SWPC, NASA OMNIWeb, ESA Swarm. The `Data/Inputs` directory currently uses mock data.
7. **Units matter**: The ESM uses real physics units (W/m², degrees C, GtCO2). Maintain unit consistency in all calculations.
