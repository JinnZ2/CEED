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

All five still execute. Run them to see the falsified behaviour directly:

```bash
cd legacy && python convergence_model_v0.py
```

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
