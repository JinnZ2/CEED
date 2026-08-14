# Legacy — the falsification record

Superseded model versions, kept because **precedence carries**.

A hypothesis that failed is a result. Deleting the code that encoded it
destroys the evidence and invites the next person — human or model — to
propose it again. This directory holds the falsified versions in runnable
form, each linked to the test that killed it and the revision that replaced
it.

The loop this project runs:

```
hypothesize -> run -> result -> falsified? -> edit the claim
                                           -> search for unknowns
                                           -> rerun
```

Everything here has been through it at least once.

## What is in here

| File | Encoded hypothesis | Falsified by | Status |
|---|---|---|---|
| `convergence_model_v0.py` | Retention as `1.1*E`; absolute phase thresholds | E1, E2 | superseded |
| `convergence_model_extended_v0.py` | Stochastic events as a separate model | E6 | merged into the main model |
| `CEED_universal_model_v0.py` | Feedback as a state multiplier; retention collapse | E5 | superseded |
| `minimum_esm_code_v0.py` | Retention collapse as runaway safeguard; ECS as forcing | E3, E4 | superseded |
| `run_mc_v0.py` | Same ECS-as-forcing error, sampled | E4 | superseded |
| *(no file)* | "Mean R² is positive" as a claim of skill | E8 | claim retracted |
| *(no file)* | "Add an ENSO oscillator" as the fix for E8 | E9 | wrong shape; replaced by `simulation/tipping.py` |
| *(no file)* | Cascade coupling to raw `h` | E10 | inverted stabilising links; fixed in `cascade.py` |
| *(no file)* | "Heavy tails always dominate" | E11 | true only when the fold is far from typical event size |

All five files still execute. Run them to see the falsified behaviour
directly:

```bash
cd legacy && python convergence_model_v0.py
```

E8 has no file because what failed was a *claim about the model*, not a
version of it. It is recorded here anyway — a retracted claim is as much a
result as a retired implementation.

## What is deliberately NOT in here

The original files at commit `6cb8779` do not parse. Four of six carried
markdown code fences, curly quotes, and markdown-bolded `**name**` from being
pasted out of a chat transcript, with indentation flattened outside the
fences. That is transcription damage, not a hypothesis — there is nothing to
falsify and nothing to learn beyond "check that pasted code compiles." It
stays in git history at `6cb8779` and does not need a runnable copy.

The distinction matters: **a broken file is not a failed experiment.**

---

## E1 — Retention as a level multiplied into a rate

**Hypothesis.** Writing retention as `1.1 * E * (1 - lambda)` means the system
keeps slightly more energy than it sheds, so accumulation emerges from the
dynamics.

**Test.** Linearise the per-subsystem derivative; integrate 3 years.

**Result — falsified.** Collecting linear coefficients:

| System | lambda | net linear coeff | implied equilibrium |
|---|---|---|---|
| solar | 0.05 | **+0.9950 /yr** | E ≈ 995 |
| magnetic | 0.02 | **+1.0580 /yr** | E ≈ 1058 |
| atmospheric | 0.08 | **+0.9320 /yr** | E ≈ 932 |
| oceanic | 0.01 | **+1.0790 /yr** | E ≈ 1079 |

Every subsystem grows at roughly **+100% per year**, checked only by
`0.001*E²`. Total energy runs 500.5 → 3021 in three years. The model cannot
express a stable outcome for *any* input.

Runaway was not a finding the model produced. It was imposed by the
coefficient. `1.1` was a retained *fraction* being used as a growth
multiplier.

**Revision.** Retention became an explicit per-subsystem `alpha`, separate
from dissipation `lambda` and quadratic loss `gamma`. Net coefficients are now
negative (solar: `alpha - lambda = -0.29 /yr`), giving real equilibria. A
10-year run goes 500.5 → 556.1.

---

## E2 — Absolute phase thresholds

**Hypothesis.** Total energy crossing 120 / 150 / 200 / 300 marks four
escalating phases of system stress.

**Test.** Classify a 3-year baseline run.

**Result — falsified.** The baseline total is
`180 + 92.5 + 118 + 110 = 500.5`, already **200.5 above the Phase 4
threshold** before a single timestep. All 36 timesteps classify as Phase 4.
Phases 1–3 are unreachable branches. `phase_1: 120` is never compared against
anything at all.

The same defect propagated: `analyze_event_timing` keyed "critical" to
`phase_2`, so a 3-year run flagged **1,756 of 1,756** events as critically
timed — informationally identical to flagging none.

**Revision.** Thresholds became multiples of the baseline total
(1.20 / 1.50 / 2.00 / 3.00 × `E_baseline`), which is what 120/150/200/300 read
as naturally — percentages. A Phase 0 band covers at-or-below baseline.
`analyze_event_timing` now measures against baseline with an explicit
`threshold_ratio`.

**Unknown surfaced.** Fixing the scale did not fix the *quantity*. The total
is still a sum of unlike indices — an F10.7 proxy added to an ocean-heat
proxy. The ratio is now meaningful; the sum still is not. Open as finding 4
in [`../Docs/numerical-audit.md`](../Docs/numerical-audit.md).

---

## E3 — Retention collapse as the runaway safeguard

**Hypothesis.** `retention = 1.05 * exp(-0.0008 * T^2)` falls below 1.0 at
high temperature, preventing unrealistic runaway. Labelled `CRITICAL FIX` in
the source, with comments claiming 0.98 at T=3 and 0.92 at T=6.

**Test.** Evaluate the function.

**Result — falsified.**

| T (°C) | comment claims | actual |
|---|---|---|
| 0 | 1.05 | 1.0500 |
| 3 | 0.98 | **1.0425** |
| 6 | 0.92 | **1.0202** |

Retention does not cross 1.0 until **T = 7.8 °C**. Across the entire range the
model realistically explores it stays above 1.0 — accumulating. The stated
safeguard never engaged.

Matching the comments would need the rate 4.6× to 9.6× larger, and those two
anchor points are **mutually inconsistent** under an `exp(-k·T²)` form — so
the intended curve could not be recovered from the comments. Either the
constant was wrong by an order of magnitude or the comments were; they could
not both be right.

**Revision.** The retention multiplier was removed entirely in favour of the
standard energy-balance form `C dT/dt = F_total - lambda_eff * T` with
`lambda_eff = F_2xCO2 / ECS = 1.23 W/(m² K)`. Bounded behaviour now comes from
the feedback parameter rather than an ad-hoc multiplier.

---

## E4 — A temperature sensitivity used as a radiative forcing

**Hypothesis.** `CO2_forcing = ECS * log2(CO2/CO2_0)`, summed with aerosol ERF
and solar forcing.

**Test.** Dimensional check at a CO₂ doubling.

**Result — falsified.** ECS has units of **°C per doubling**, not W/m².

| | value at 2×CO₂ |
|---|---|
| as coded, `ECS*log2(C/C0)` | **3.00** (this is ECS, in °C) |
| standard, `5.35*ln(C/C0)` | **3.71 W/m²** |
| AR6 assessed | ~3.93 W/m² |

A temperature was being summed with three fluxes and fed through retention and
dissipation as though it were an energy flux. The error was masked by
coincidence: ECS ≈ 3 and forcing-per-doubling ≈ 3.7–3.9 are close enough that
output looked plausible.

**Revision.** `F_2xCO2 = 3.7 W/m²` (Myhre et al. 1998) is now the forcing
constant; ECS enters only through `lambda_eff = F_2xCO2 / ECS`. Applied in
both the ESM and the Monte Carlo sampler.

---

## E5 — Feedback as a state multiplier

**Hypothesis.** Each feedback rewrites the state, `E -> E*(1 ± strength)`,
with retention `1.05 * exp(-0.001 * E^2)` and thresholds at 100 / 200 / 300.

**Test.** Tabulate retention against the thresholds; drive with escalating
forcing.

**Result — falsified.** Retention is effectively zero long before the
*warning* threshold:

| E | retention |
|---|---|
| 50 | 0.0862 |
| 100 | **4.77e-05** |
| 300 | **8.60e-40** |

The shipped demo peaks at energy **51.9** against a warning threshold of 100.
Buffer depletion is gated on `energy > warning_threshold`, so
`buffer_capacity` stayed pinned at 100% and `stressed` / `critical` /
`tipping` were unreachable. A forcing of 5000 peaked at 1551 and still
reported `stable`.

Also incoherent internally: retention and dissipation were computed from
`self.state.energy` while feedbacks were applied to a separate running
`net_effect`, so the two halves of `update()` operated on different
quantities.

**Revision.** Feedbacks became rate contributions with rational saturation,
`f_k = ±s_k·E/(1 + (E/E_sat)²)`, and the system integrates
`dE/dt = F_ext + Σ f_k − D(E)`. Thresholds are now reachable: a constant
forcing of 50 drives the system to peak energy 474 and classifies as
`tipping`.

---

## E6 — Stochastic calls inside the ODE right-hand side

**Hypothesis.** External events and uncharacterised sinks belong in a separate
extended model.

**Test.** Inspect the derivative function.

**Result — falsified on correctness, not on physics.** `unknown_dissipation`
called `np.random.random()` **inside** the ODE right-hand side and appended to
`self.events` during integration. `odeint` evaluates the derivative several
times per step and assumes it is a deterministic function of `(E, t)`; a
stochastic RHS silently corrupts the solver's error control and step-size
selection. The model also duplicated the E1 retention error.

**Revision.** Folded into `ConvergencePredictor` as a subclass in the same
file, with all events **pre-generated before integration** so the RHS stays
deterministic. Guarded by
`test_events_are_pregenerated_so_the_rhs_is_deterministic`.

---

## E7 — Solar cycle phased on the source instead of the state

The cleanest complete turn of the loop, and the only entry whose falsification
came from the hindcast rather than from inspection.

**Hypothesis.** Phase the solar source to peak at the observed F10.7 maximum
(2014.5) and the model will reproduce the 11-year cycle. Committed as
*"Recalibrate model against 2010-2024 hindcast: mean R² now positive."*

**Test.** Hindcast against F10.7 annual means, 2010–2024.

**Result — falsified.** **R² = −0.443, grade F** — worse than predicting a
flat mean. Predicted peak 2016–2017 against an observed peak of 2014, with
errors reaching −74 sfu. The reported positive *mean* R² (+0.183) was carried
entirely by two other subsystems and hid this one.

**Search for unknowns.** F10.7 is the *state* `E_solar`, not the source, and
the subsystem is a first-order low-pass filter. Linearising about equilibrium:

```
lambda_eff = (lambda - alpha) + 2*gamma*E_eq  = 0.470 /yr
lag        = arctan(omega/lambda_eff)/omega   = 1.545 yr
gain       = 1/sqrt(1 + (omega/lambda_eff)^2) = 0.635
```

Phasing the source at the observed maximum necessarily put the state peak
1.5 years late. The same filter attenuates amplitude by 0.635, so the
comment's claim that ±0.30 modulation reproduced the 70–180 sfu range was also
wrong.

**Edited claim.** Source must *lead* by the lag: offset `0.5 → 2.045`.
Modulation pre-compensated for the gain: `0.30 → 0.485`. Both derived from the
model's own linearisation — **not fitted to the hindcast**, so the hindcast
remains an independent check.

**Rerun.** **R² = +0.359, grade C.**

---

## E8 — "Mean R² is positive" as a claim of skill

No superseded *file* for this one. What was falsified is a **claim about the
model**, made repeatedly in commit messages and docs, including by this
project's own audit.

**Hypothesis.** The model has predictive skill, evidenced by a positive mean
R² against 2010–2024 observations — first +0.183, then +0.312 after the E7
correction.

**Test.** Split the record. Train 2010–2019, test 2020–2024, one continuous
trajectory with no re-initialisation at 2020, so the test years are a
free-running five-year forecast. Re-derive the solar amplitude from the train
window alone, since the full-record derivation reads the 2024 F10.7 maximum
and that sits in the test window.

**Result — falsified.** Shipped parameters, held out:

| System | train R² | test R² | gap |
|---|---|---|---|
| solar | −0.020 | **+0.555** | −0.575 |
| magnetic | −1.193 | +0.057 | −1.249 |
| atmospheric | +0.316 | +0.042 | +0.274 |
| oceanic | **+0.773** | **−2.705** | **+3.478** |
| **mean** | −0.031 | **−0.513** | +0.482 |

Mean test R² of **−0.513**: out of sample the model is beaten by a flat line
through the test window's own mean. The reported +0.312 was measuring the
optimiser, not the physics.

**Oceanic is the specific failure.** The subsystem that looked strongest
in-sample (+0.631, grade B) fails hardest out of sample, a gap of 3.478. It
was fitting a monotonic trend it cannot extrapolate. Refitting without
leakage halves the damage (−0.469) but does not remove it.

**Solar is the counterexample worth keeping.** It is the only subsystem with a
*negative* gap under both parameter sets — better out of sample than in. Its
constants came from the model's own linearisation (E7), not from fitting.
Derived constants travelled; fitted ones did not. That is the most useful
thing this test produced.

**Edited claim.** Full-record scores are labelled in-sample everywhere they
appear, and the held-out table is now the one quoted as skill. The
calibration guide's old rule — *"commit if scores improved"* — is itself
falsified by this result and replaced: commit only if the **test** score
improved, or if the change is derived rather than dialled.

**Unknown surfaced.** The test window is 5 years against an 11-year cycle, so
it cannot test cycle physics at all — only drift, bias, and extrapolation. A
real test of the coupling matrix needs either a longer record or
cross-system lag correlations, which nothing here scores yet.

Reproduce with `python -m simulation.validation`. The separation is verified,
not assumed: `test_calibration_cannot_see_the_test_window` poisons every
test-window observation and asserts the fitted parameters do not move.

---

## E9 — "Add an oscillator" as the fix for the oceanic failure

Also no file. This entry records a hypothesis *this project proposed and then
killed within the same session*, which is the cheapest kind of falsification
and worth keeping visible.

**Hypothesis.** The oceanic subsystem's out-of-sample collapse (E8) is caused
by a missing ENSO mode. Adding a recharge oscillator — the Jin (1997)
two-equation form — would supply the interannual variability the model cannot
produce, and the fit would improve.

**Supporting evidence, which was real.** Model residuals are not white:
lag-1 autocorrelation +0.917 oceanic, +0.839 solar, +0.363 atmospheric.
Detrended observed temperature matches documented ENSO phase in **9 of 11**
non-neutral years, with the largest excursions exactly where expected (2016
peak +0.161, the 2021–22 La Niña −0.177/−0.172, 2024 peak +0.147). The
detrended ENSO range is 0.338 K against a model atmospheric RMSE of 0.13 K —
the signal the model cannot make is 2.6× its total error.

**Falsified — wrong shape, not wrong topic.** A recharge oscillator produces
variance around a *fixed* equilibrium. It cannot produce a state that does not
come back. But the observed record shows exactly that: Antarctic sea ice broke
in September 2016 and the 2023–2025 record lows sit inside the new regime, not
as excursions from the old one. [Nature Communications
(2025)](https://www.nature.com/articles/s41467-025-66143-7) finds super El
Niños drive "abrupt, persistent transitions ... for years or even decades"
after the event fades. An oscillator models none of that.

**Second sub-hypothesis, also falsified.** ENSO is an ocean↔atmosphere heat
exchange, so the two should be *anti*-correlated — the discharge signature, and
a direct test of CEED's coupling term. Measured: `corr(T, OHC) = +0.338`, with
El Niño years showing *higher* upper-ocean heat content. Not supported. Global
0–700 m OHC is trend-dominated, the discharge signal is tropical-Pacific
specific, and annual resolution smears ENSO's ~6-month lead/lag. This repo's
OHC series is also coarse — 0.5×10²² J steps, with a 2017–19 dip that looks
like an artifact.

**Edited claim.** The requirement is not an oscillator but a **fast–slow
bistable element**: two branches, a fold, hysteresis, and separated
timescales so commitment precedes visible change. Built as
`simulation/tipping.py`, and deliberately kept standalone — the convergence
model still has one attractor, and wiring the two together changes the model's
character rather than fixing a defect.

**Unknown surfaced.** The regime break is September 2016, which sits *inside*
the 2010–2019 training window. That is a better explanation of E8's oceanic
failure than a missing oscillator: a stationary model fitted across a
non-stationary record produces a decent in-sample blend and a bad
extrapolation. Fitting pre- and post-break separately does shift the
parameters (`A_0` 2.94 → 4.78), but `atm_fraction` hit the search boundary on
~5 scored points, so this repo's 15 annual values **cannot** establish the
break. The sea-ice literature establishes it with 45 years of daily data; take
it from there, not from here.

---

## E10 — Coupling a cascade to the raw slow variable

A bug caught by its own test within minutes of being written, kept because
the failure mode is subtle and would have looked like a modelling result.

**Hypothesis.** In a network of coupled tipping elements, element j's
influence on element i is `C_ij * h_j`, where h is j's slow observable.

**Test.** Give the link a negative (stabilising) coefficient and check that
the target does not tip.

**Result — falsified.** The target tipped at t=2.0, almost immediately. `h`
runs from −1 (intact) to +1 (gone), so at t=0 every element has h = −1 and a
negative coefficient times −1 is a **positive** shove. Two bugs in one:

1. An untouched network exerted influence before anything had tipped.
2. Stabilising links were inverted — they *caused* the cascade they were
   meant to prevent.

The sign error is the dangerous one. It would not have crashed anything. It
would have produced a plausible-looking figure showing that protective
couplings accelerate collapse, which is the opposite of what the model says.

**Revision.** Couple to the tipped fraction `phi_j = (h_j + 1) / 2`, which is
0 when intact and 1 when fully tipped. Influence is now zero until something
actually tips, and a negative coefficient is protective. Guarded by
`test_untouched_network_exerts_no_influence` and
`test_stabilising_link_is_protective_not_inverted`.

**Unknown surfaced.** Fixing it made the original demo of a protective link
stop working — a stabilising Ross→Thwaites coupling of *any* strength failed
to save Thwaites. That is not a second bug. Ross needs `tau_slow = 400` to
respond and Thwaites tips at `t = 4`, so the help arrives four hundred years
late. The demo now shows that instead, because it is the more useful result:
**a protective coupling slower than the collapse it prevents is worthless.**
Sign alone does not determine whether an interaction helps.

---

## E11 — "Heavy tails always dominate the tipping risk"

**Hypothesis.** Rare extreme events tip a bistable system more readily than
ordinary fluctuations of the same mean and variance, because the tail reaches
further. Predicted before running: the heavy-tailed case would dominate at
every forcing level.

**Test.** Poisson event trains with identical mean and identical variance —
Exponential(mu) against Normal(mu, mu) — driving a tipping element whose mean
forcing is held below the fold.

**Result — falsified as stated.** At mean magnitude 0.30 the ordering was
*backwards*: thin tail 89%, heavy tail 75% at F=0.

**Diagnosis.** The fold sat only 1.28x the typical event size away. At that
distance Normal(mu, mu) puts *more* mass past the threshold (0.389) than
Exponential(mu) does (0.277), because the exponential piles probability near
zero to pay for its tail. The tail only wins further out.

| mean magnitude | fold / mean | P(exp > fold) | P(norm > fold) | ratio |
|---|---|---|---|---|
| 0.30 | 1.28 | 0.2772 | 0.3886 | **0.7x** |
| 0.20 | 1.92 | 0.1459 | 0.1776 | 0.8x |
| 0.15 | 2.57 | 0.0768 | 0.0587 | 1.3x |
| 0.10 | 3.85 | 0.0213 | 0.0022 | **9.7x** |
| 0.06 | 6.42 | 0.0016 | 0.0000 | 53000x |

**Edited claim.** Heavy tails dominate *when the threshold is far out relative
to the typical event*. Re-run at mean magnitude 0.10 (fold 3.85x away), the
predicted ordering appears: 27% heavy against 15% thin at F=0.

**What survived, and it is the more important half.** The headline result was
never about tails. At F=0 — barrier at its maximum, mean forcing nowhere near
the fold — the system tips in 15-27% of 400-year trials regardless of
distribution. **A threshold is not a safety margin once discrete events are in
the picture.** That claim was strengthened, not weakened, by the tail
hypothesis failing.

**Unknown surfaced.** Reporting event risk as a variance understates it when
the threshold is far out and overstates it when close. Neither the convergence
model nor the red team report distinguishes these regimes; both speak of
"thresholds" as though crossing required the mean to arrive.

---

## Falsified but not yet revised

Open claims, recorded here so they are not mistaken for settled. Full detail
in [`../Docs/numerical-audit.md`](../Docs/numerical-audit.md).

| Claim | Where | Status |
|---|---|---|
| "Dipole field weakening 9–10% per decade" | `red-team-report.md` | Measured rate is ~5% per century; 9% is the decline since **1840**. Off by ~10×. Correction pending. |
| "Thermosphere still expanded from 2024 storm" | `red-team-report.md` | Storm density recovers over days, not 14 months. Correction pending. |
| "AMOC slowdown faster than models suggest" | `red-team-report.md` | RAPID shows −0.8 ± 0.7 Sv/decade — uncertainty nearly as large as the trend. Overstated. |
| Scenario probabilities (40% / 25% / 30% / 10% / 5%) | `red-team-report.md` | No horizon, method, or source. |
| Magnetic subsystem tracks the solar cycle | `convergence_model.py` | **R² = −0.214**, grade F. Source is constant; solar coupling carries only ~2% of its input, so it cannot oscillate. Structural, not a tuning error. |

The magnetic entry is deliberately unfixed. Tuning it against the same
2010–2024 window used to score it would produce a better number and no new
knowledge.

---

## Provenance

| Commit | What it was |
|---|---|
| `6cb8779` | Original state; four of six Python files did not parse |
| `e2606b4` | Transcription damage repaired, original designs intact — **source of the v0 files here** |
| `df72ed6` | MHD/spatial module and first numerical audit |
| `2f1ce69` | Independent model rewrite (hindcast, unit bridge, alpha retention) |
| `16bd52f` | The two lines merged; E7, E2 and the normalisation defect fixed |

Retrieve any file at any point with:

```bash
git show <commit>:<path>
```
