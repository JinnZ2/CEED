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
python -m simulation.validation             # held-out train/test split
```

## Current Scores

**In-sample** (full record 2010-2024, the window parameters were tuned on):

| System      | Baseline R² | Anthro R² | Grade |
|-------------|-------------|-----------|-------|
| Solar       | -0.45       | +0.36     | C     |
| Magnetic    | -0.02       | -0.21     | F     |
| Atmospheric | -0.22       | +0.47     | C     |
| Oceanic     | -1.45       | +0.63     | B     |
| **Mean**    | -0.53       | **+0.31** |       |

**Out-of-sample** (train 2010-2019, test 2020-2024) — this is the honest one:

| System      | train R² | test R²    | gap    |
|-------------|----------|------------|--------|
| Solar       | -0.020   | **+0.555** | -0.575 |
| Magnetic    | -1.193   | +0.057     | -1.249 |
| Atmospheric | +0.316   | +0.042     | +0.274 |
| Oceanic     | +0.773   | **-2.705** | +3.478 |
| **Mean**    | -0.031   | **-0.513** | +0.482 |

The in-sample mean of +0.31 becomes **-0.513** held out. Quote the second
table, not the first.

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

### Solar (was F, now C — and it generalises)
**Fixed.** The cycle was phased on the *source* while F10.7 is the *state*.
The subsystem is a first-order low-pass, so the state lags the source by
`arctan(omega/lambda_eff)/omega` = 1.545 yr and is attenuated by 0.635.
`solar_phase_offset` and `solar_modulation` now carry the compensation, both
derived rather than fitted.

This is the only subsystem with a **negative** train-test gap, i.e. the only
one that does better out of sample than in. That is what a derived constant
buys you.

**Remaining**: a single cosine still cannot capture the asymmetric rise/fall
of real cycles, or the differing amplitudes of SC24 and SC25. Harmonics or a
data-driven source would help — but note that fitting harmonics to 15 years
of data is exactly the move the held-out test is there to catch.

### Magnetic (F grade)
**Problem**: the model produces a smooth monotonic Kp (1.30 -> 1.91) while
observations oscillate 1.3-2.3. The source is a constant 17.7 energy/yr, so
all variability must arrive through the solar coupling — which delivers
`c*eta*E_solar = 0.015 x 0.15 x 180 ~ 0.4` energy/yr against that source,
about **2%**. It is structurally unable to oscillate.

**Next step**: a solar-modulated magnetic source term, or much stronger
solar->magnetic coupling. Justify either against data — do not fit it to
2010-2024 and then report a 2010-2024 score.

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
# 3. Run the in-sample hindcast
python -m simulation.hindcast --compare
# 4. Run the held-out split — THIS is the one that decides
python -m simulation.validation
# 5. If atmospheric R² dropped, you broke the energy budget
# 6. If oceanic R² dropped, coupling or A_0 is off
# 7. Commit only if the TEST score improved, or if the change is
#    derived from physics rather than fitted
```

**"Commit if scores improved" was the old rule and it was wrong.** Improving
the in-sample score is trivial — add parameters. The oceanic subsystem scores
+0.773 in-sample and **-2.705** held out, which is what that rule produces if
you follow it long enough.

Two questions before committing a parameter change:

1. Did the **test** score improve, or only the train score?
2. Is the new value **derived** from something (a timescale, a conservation
   law, a published measurement) or was it dialled until the number looked
   good? Derived constants generalise. The solar fix is the worked example:
   derived from the model's own linearisation, and the only subsystem with a
   negative train-test gap.

A parameter fitted on 2010-2024 and then scored on 2010-2024 tells you
nothing except that the optimiser worked.
