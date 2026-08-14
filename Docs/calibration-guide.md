# CEED Calibration Guide

How to tune the convergence model against observations.

## What This Model Is

CEED is **not** a replacement for full Earth System Models (CESM, GFDL, UKESM).
Those models solve the primitive equations on 3D grids with resolved physics.
CEED can't compete on single-system accuracy.

**CEED's thesis**: Earth systems (solar, magnetic, atmospheric, oceanic) are
coupled through energy transfer pathways that existing models treat in
isolation.  The model's value is in the **coupling matrix** — testing whether
cross-domain energy transfer can amplify anthropogenic forcing beyond what
single-system models predict.

The hindcast is the ground truth.  If the model can't reproduce the past,
it can't predict the future.

## Running the Hindcast

```bash
python -m simulation.hindcast --compare     # baseline vs anthropogenic
python -m simulation.hindcast --system solar  # single system detail
python -m simulation.hindcast --plot        # with plots
```

## Current Scores (as of calibration v2)

| System      | Baseline R² | Anthro R² | Grade | Status              |
|-------------|-------------|-----------|-------|---------------------|
| Solar       | -0.45       | -0.44     | F     | Needs work          |
| Magnetic    | -0.02       | -0.03     | F     | Near zero (neutral) |
| Atmospheric | -0.22       | +0.51     | C     | Calibrated          |
| Oceanic     | -1.45       | +0.69     | B     | Calibrated          |
| **Mean**    | -0.53       | **+0.18** |       | Positive            |

## Parameter Hierarchy

Tune parameters in this order.  Each level depends on the previous being
correct.

### Level 1: Equilibrium (source rates)

Each subsystem's source rate sustains its equilibrium energy.  At equilibrium:

```
S = (lambda - alpha) * E_eq + gamma * E_eq^2
```

| System      | E_eq  | S_needed | Current S | Status |
|-------------|-------|----------|-----------|--------|
| Solar       | 180.0 | 68.4     | 68.4      | Set    |
| Magnetic    | 92.5  | 17.7     | 17.7      | Set    |
| Atmospheric | 118.0 | 2.55     | 2.55      | Set    |
| Oceanic     | 110.0 | 2.20     | 2.20      | Set    |

**Rule**: if you change alpha, lambda, or gamma, recalculate S to maintain
the same equilibrium.  Otherwise the system drifts.

### Level 2: Response timescale (lambda)

Lambda controls how fast each subsystem responds to perturbations:

| System      | lambda | Response time (1/λ) | Physical basis            |
|-------------|--------|---------------------|---------------------------|
| Solar       | 0.30   | ~3 years            | F10.7 tracks current activity |
| Magnetic    | 0.15   | ~7 years            | Kp responds to solar wind |
| Atmospheric | 0.08   | ~12 years           | Tropospheric thermal inertia |
| Oceanic     | 0.01   | ~100 years          | Deep ocean mixing time    |

**Constraint**: solar lambda must be large enough that the 11-year cycle
isn't damped out.  Currently ±68 sfu amplitude (matching observations)
requires lambda ≥ 0.2.

### Level 3: Anthropogenic forcing (A_0, atm_fraction)

Controls how much geological energy release goes to each system:

| Parameter    | Current | Constraint                                  |
|--------------|---------|---------------------------------------------|
| A_0          | 2.5     | Total energy injection rate                 |
| atm_fraction | 0.22    | ~22% to atmosphere, ~78% to ocean           |
|              |         | Ocean absorbs ~93% of excess heat (AR6)     |

**Tuning**: scan A_0 from 0.5 to 5.0 and atm_fraction from 0.15 to 0.50.
Score atmospheric and oceanic R² jointly.  There is a Pareto tradeoff:
higher A_0 helps oceanic but hurts atmospheric (and vice versa).

### Level 4: Coupling matrix

9 pathways with (rate, efficiency) pairs.  Tune only after levels 1-3
are set.

**Priority couplings** (most impact on hindcast):
1. `(magnetic, solar)`: controls whether Kp tracks solar cycle
2. `(oceanic, atmospheric)`: controls ocean heat uptake rate
3. `(atmospheric, oceanic)`: controls ENSO-like atmospheric feedback

**Constraint**: all efficiencies must be in (0, 1).  Coupling-only system
must lose total energy (waste heat test).

### Level 5: State-dependent parameters

Only tune after static parameters are working:
- `alpha_sensitivity`: how much atmospheric retention increases with CO2
- `coupling_sensitivity`: how gradients strengthen coupling
- `efficiency_shift`: how total energy shifts conversion efficiency

## Known Issues and Next Steps

### Solar (F grade)
**Problem**: the model tracks the cycle direction but amplitudes don't match.
The cos-based source function captures the period but not the asymmetric
rise/fall of real solar cycles.  SC24 and SC25 have different amplitudes.

**Next step**: replace single cosine with superposition of harmonics, or
use a data-driven source function that takes observed F10.7 as input for
hindcast and a parametric cycle for projection.

### Magnetic (F grade, near zero)
**Problem**: R²≈0 means the model produces roughly constant Kp while
observations show solar-cycle-correlated variability.  The solar→magnetic
coupling (c=0.015, η=0.15) transfers energy too slowly relative to
magnetic dissipation (λ=0.15).

**Next step**: strengthen solar→magnetic coupling or reduce magnetic
dissipation so that solar cycle variability propagates into Kp.

### Cross-system validation
**Problem**: current scoring treats each system independently.  CEED's
thesis is about cross-system correlations.

**Next step**: add lag-correlation scores between subsystem pairs
(e.g., does F10.7 peak predict Kp peak with the right lag?).  This
would test the coupling matrix directly.

## Quick Calibration Workflow

```bash
# 1. Change a parameter in convergence_model.py
# 2. Run tests
python -m pytest -q
# 3. Run hindcast
python -m simulation.hindcast --compare
# 4. Check scores — did they improve?
# 5. If atmospheric R² dropped, you broke the energy budget
# 6. If oceanic R² dropped, coupling or A_0 is off
# 7. Commit if scores improved
```
