# CEED Numerical Audit

**Scope:** parameters, thresholds, units, and literature-facing values across
the model code and the risk documentation.
**Method:** every internal finding below was computed directly from the
shipped code, not estimated. Literature values were checked against primary
or near-primary sources, cited inline.
**Status:** revised after the model-rewrite branch was merged. Findings that
the rewrite or this merge resolved are marked **RESOLVED** with what fixed
them; the rest are open. Reproduce any figure with `pytest` and the module
`__main__` blocks.

---

## Summary

| # | Finding | Where | Severity | Status |
|---|---------|-------|----------|--------|
| 1 | Phase thresholds unreachable; every run Phase 4 from t=0 | `convergence_model.py` | Critical | **RESOLVED** |
| 2 | Retention term made runaway structural | `convergence_model.py` | Critical | **RESOLVED** |
| 3 | Retention-collapse "fix" never engaged | `minimum_esm_code.py` | Critical | **RESOLVED** |
| 4 | Incommensurable units summed into one total | `convergence_model.py` | High | **Partly resolved** |
| 5 | ECS (°C/doubling) used where a forcing (W/m²) was required | `minimum_esm_code.py` | High | **RESOLVED** |
| 6 | Universal model thresholds outside reachable range | `CEED_universal_model.py` | High | **RESOLVED** |
| 7 | Dipole decay rate overstated by ~10x | `red-team-report.md` | High | **Open** |
| 8 | Roadmap formulas not dimensionally closed | `To-be-added.md` | Medium | Documented |
| 9 | Declared parameters never read | `minimum_esm_code.py` | Medium | **RESOLVED** |
| 10 | Emissions rate contradicted its own comment | `minimum_esm_code.py` | Low | **RESOLVED** |
| 11 | Sink split drifted from Global Carbon Budget | `minimum_esm_code.py` | Low | **RESOLVED** |
| 12 | Unsourced probabilities; two overstated observations | `red-team-report.md` | Medium | **Open** |
| 13 | Solar cycle phased on the source, not the state | `convergence_model.py` | High | **RESOLVED** |
| 14 | Conversion-efficiency response missing its normalisation | `convergence_model.py` | Medium | **RESOLVED** |
| 15 | Hindcast scored the year it was initialised from | `hindcast.py` | Medium | **RESOLVED** |
| 16 | Ap→Kp conversion did not reproduce its own table | `unit_bridge.py` | Medium | **RESOLVED** |
| 17 | `pytest` collected zero tests | `Tests/` | High | **RESOLVED** |
| 18 | Magnetic subsystem cannot track the solar cycle | `convergence_model.py` | Medium | **Open** |
| 19 | Reported skill was in-sample; held-out test is far worse | `hindcast.py` | **Critical** | **Open** |

---

## Resolved in the model rewrite

**1 — Phase thresholds.** Thresholds are now multiples of the baseline total
(`phase_threshold_ratios`, 1.20/1.50/2.00/3.00 × `E_baseline` = 500.5) rather
than absolute values of 120/150/200/300 that sat below the baseline. A new
Phase 0 band covers "at or below baseline". A 3-year baseline run now
classifies as Phase 0 rather than Phase 4, and
`test_baseline_run_is_not_pinned_at_phase_four` guards the regression.

**2 — Structural runaway.** Retention is now an explicit `alpha` rate per
subsystem, separate from dissipation `lambda` and quadratic loss `gamma`,
rather than the old `1.1*E*(1-lambda)` level-as-rate term. Net linear
coefficients are now negative (e.g. solar: `alpha - lambda = -0.29 /yr`), so
the model has genuine equilibria. A 10-year run goes 500.5 → 556.1 instead
of 500.5 → 3021.

**3, 5, 9 — Minimal ESM.** The model was rebuilt on a standard energy-balance
form: `C dT/dt = F_total - lambda_eff*T` with `lambda_eff = F_2xCO2/ECS =
1.23 W/(m² K)`. The broken `retention_factor` (which never dropped below 1.0
until 7.8 °C, contradicting its own comments) is gone. CO₂ forcing is now
`F_2xCO2 * log2(CO2/CO2_0)` with `F_2xCO2 = 3.7 W/m²` (Myhre et al. 1998)
instead of ECS standing in for a flux. All 15 declared parameters are read;
previously 6 were dead.

**6 — Universal model.** Feedbacks now use a saturating form
`f_k = ±s_k·E/(1+(E/E_sat)²)` and thresholds are reachable: a constant
forcing of 50 drives the system to a peak energy of 474 and a `tipping`
classification. Previously a forcing of 5000 peaked at 1551 and still
reported `stable`, so buffer depletion never ran.

**10, 11 — Carbon numbers.** `emissions_rate` is 11.1 GtC/yr, matching
[Global Carbon Budget 2024](https://essd.copernicus.org/articles/17/965/2025/)
(11.1 GtC/yr in 2023, 11.3 preliminary for 2024); it was 10.0 under a comment
that read "current ~11 GtC/yr". The incorrect land/ocean sink split (0.31 /
0.23 against an observed 0.29 / 0.26) is gone — only the total, 0.54, remains,
and it is used.

---

## Resolved in this merge

**13 — Solar cycle phased on the wrong quantity.** This was the single largest
error in the hindcast. F10.7 is the *state* `E_solar`, not the source, and the
subsystem is a first-order low-pass filter. Linearising about equilibrium:

```
lambda_eff = (lambda - alpha) + 2*gamma*E_eq = 0.470 /yr
lag        = arctan(omega/lambda_eff)/omega  = 1.545 yr
gain       = 1/sqrt(1 + (omega/lambda_eff)^2) = 0.635
```

The source was phased to peak at the observed 2014.5 F10.7 maximum, which put
the *state* peak at 2016.0. The source must lead by the lag, so the offset
goes from 0.5 to **2.045**. The same filter attenuates amplitude by 0.635, so
the modulation goes from 0.30 to **0.485** to reproduce the observed ±55.5 sfu
half-amplitude. Both constants are derived, not fitted.

Effect on the hindcast: solar **R² −0.443 (F) → +0.359 (C)**.

**14 — Missing normalisation.** `effective_coupling` computed
`eta_eff = eta_base*(1 + delta*(E_total - E_ref))`, omitting the `/E_ref` its
own docstring specified. That inflated the response by a factor of
`E_ref` = 500.5 and left `delta` carrying units of 1/energy. The division is
restored and `efficiency_shift` rescaled 0.0003 → 0.15 (= 0.0003 × 500.5) so
calibrated behaviour is unchanged while the expression is dimensionless.

**15 — Self-flattering hindcast.** The model is initialised from the observed
value at `start_year`, so its residual there was identically zero. That year
was also being scored, donating one free perfect point per subsystem. Scoring
now starts at `start_year + 1`.

**16 — Ap→Kp conversion.** `unit_bridge.py` claimed `Kp ≈ 0.3*Ap^0.55`
generated its Kp table. It does not — it is low by up to 1.07 Kp units across
the record, and also falls below the official GFZ ap↔Kp scale (ap=5 → Kp 1.50,
ap=9 → 2.25, ap=14 → 2.88). The claim is removed; the table is documented as
published annual means of 3-hourly Kp, with the caveat that averaging a
quasi-logarithmic index is itself a compromise. `OBSERVED_AP` is retained for
anyone redoing it properly.

**17 — Zero-test collection.** 22 tests existed but lived in
`Tests/test-convergence-model.py`; pytest's default discovery matches
`test_*.py`, so a bare `pytest` reported "no tests ran" and passed. Renamed,
plus `conftest.py` and package markers so `simulation.*` imports resolve and
`python -m simulation.hindcast` is no longer the only way to run it. The suite
now collects 141 tests.

Also cleaned: `years_to_model_time` (never called), the unused `OBSERVED_ENERGY`
import in `hindcast.py`, and the unused `t_model` parameter of
`energy_to_observations`.

---

## Open

### 4 — The total is still a sum of unlike indices (partly resolved)

`unit_bridge.py` is a real improvement: each subsystem now has a documented
observable, reference value, and scale factor, and conversions round-trip.
The `Kp = 3` mislabel on a value of 92.5 is gone.

What remains is that `classify_phases` still sums the four indices. Even
normalised, adding an F10.7 proxy to an ocean-heat-content proxy produces a
total whose physical meaning is undefined — the sum is dominated by whichever
subsystem happens to have the largest numeric scale. Thresholds relative to
baseline (finding 1) make the ratio meaningful without making the sum
meaningful. A weighted or explicitly index-based total would close this.

### 7 — Dipole decay overstated by ~10x

`red-team-report.md` still states:

> Earth's dipole field weakening 9–10% per decade in select zones

The measured global rate is about **5% per century**, and the axial dipole is
about **9% weaker than in 1840** — 9% over ~175 years, not per decade
([Olson & Amit](https://pages.jh.edu/polson1/pdfs/ChangesinEarthsDipole.pdf);
[Finlay et al.](https://www.nature.com/articles/ncomms10422)). The figure
matches the total decline since 1840 almost exactly, so this reads as a units
slip. The SAA does weaken faster regionally, but it is also
[not evidence of an imminent reversal](https://www.pnas.org/doi/10.1073/pnas.1722110115).

### 12 — Unsourced probabilities and two overstated observations

The scenario table (GIC Cascade 40%, AMOC Stall 25%, Nonlinear Coupling 30%,
Atmospheric Arc 10%, Crustal Shear Slip 5%) has no horizon, method, or source.

"Thermosphere still expanded from 2024 storm" (report dated July 2025):
storm-time density can rise by more than an order of magnitude but recovers
over days, not 14 months ([Parker & Linares 2024](https://arxiv.org/abs/2406.08617)).

"AMOC slowdown occurring faster than decadal models suggest": RAPID shows
**−0.8 ± 0.7 Sv/decade** since 2004 — uncertainty nearly as large as the
trend, and no clear slowdown further north
([Carbon Brief](https://interactive.carbonbrief.org/amoc-explainer/index.html)).

### 18 — The magnetic subsystem cannot track the solar cycle

Magnetic remains the one failing subsystem: **R² = −0.214**, grade F. This is
structural rather than a tuning error. Its source is a constant 17.7
energy/yr, so all cycle variability must arrive through the solar coupling —
but that path contributes only
`c·eta·E_solar = 0.015 × 0.15 × 180 ≈ 0.4` energy/yr against a source of 17.7,
about 2%. The model therefore produces a smooth monotonic rise (Kp 1.30 →
1.91) while observations oscillate between 1.3 and 2.3.

Closing this needs either a solar-modulated magnetic source term or a much
stronger solar→magnetic coupling, justified against data rather than fitted.
It was left alone here deliberately: tuning it against the same 2010–2024
window used to score it would not be evidence of anything.

### 19 — The reported skill was in-sample, and it does not hold up

Every score this project has quoted was measured on the window the parameters
were tuned against. `simulation/validation.py` splits the record — train
2010–2019, test 2020–2024 — and the result changes the story.

The split is enforced two ways. The model is initialised at 2010 and
integrated straight through 2024 without re-initialisation, so the test years
are a free-running five-year forecast rather than a warm start. And the solar
modulation, which is *derived* from the observed F10.7 amplitude, is
re-derived from the train window alone — over the full record that amplitude
includes the 2024 maximum of 180 sfu, which sits in the test window, so the
shipped constant leaks (0.485 full-record vs **0.337** train-only).

**Shipped parameters, scored out-of-sample:**

| System | train R² (2011–2019) | test R² (2020–2024) | gap |
|---|---|---|---|
| solar | −0.020 | **+0.555** | −0.575 |
| magnetic | −1.193 | +0.057 | −1.249 |
| atmospheric | +0.316 | +0.042 | +0.274 |
| oceanic | **+0.773** | **−2.705** | **+3.478** |
| **mean** | −0.031 | **−0.513** | +0.482 |

**Refitted on 2010–2019 only, no leakage:**

| System | train R² | test R² | gap |
|---|---|---|---|
| solar | +0.197 | +0.493 | −0.296 |
| magnetic | −1.191 | +0.060 | −1.250 |
| atmospheric | +0.545 | +0.137 | +0.408 |
| oceanic | +0.915 | **−0.469** | +1.384 |
| **mean** | +0.117 | **+0.055** | +0.061 |

Three things follow.

**The +0.312 full-record mean is not a skill estimate.** Held out, the shipped
parameters score a mean test R² of **−0.513** — worse, on average, than a flat
line through the test window's own mean.

**Oceanic is the overfit.** It is the subsystem that looked strongest
(+0.631, grade B on the full record) and it fails hardest out of sample:
+0.773 → **−2.705**, a gap of 3.478. Refitting without leakage roughly halves
the damage (−0.469) but does not remove it. The ocean term is fitting a
monotonic trend and cannot extrapolate it.

**Solar generalises, which is the one real result.** Its correction came from
the model's own linearisation rather than from fitting (finding 13), and it is
the only subsystem with a *negative* gap under both parameter sets — test
+0.555 shipped, +0.493 refitted. Derived constants travelled; fitted ones did
not.

Caveat on the caveat: the test window is 5 years against an 11-year solar
cycle, so this cannot test cycle physics. It tests drift, bias, and
extrapolation. `test_calibration_cannot_see_the_test_window` poisons every
test-window observation and asserts the fitted parameters do not move, so the
separation is verified rather than assumed.

Run it:

```bash
python -m simulation.validation
```

---

## Values that check out — leave these alone

| Parameter | Value | Source |
|---|---|---|
| `ECS` | 3.0 (likely 2.5–4.0) | AR6 best estimate 3 °C |
| `aerosol_ERF` | −1.1 (−1.7 to −0.4) | AR6 assessed, 5–95% range |
| `F_2xCO2` | 3.7 W/m² | Myhre et al. 1998; AR6 assesses ~3.93 — defensible, slightly dated |
| `cloud_feedback_strength` | 0.45 W/(m² K) | AR6 net cloud feedback |
| permafrost range | 14–175 GtCO₂/°C | [AR6 WG1 FAQ Ch.5](https://www.ipcc.ch/report/ar6/wg1/downloads/faqs/IPCC_AR6_WGI_FAQ_Chapter_05.pdf) |
| `total_sink_fraction` | 0.54 | GCB ~55% |
| `GtC_per_ppm` | 2.12 | standard |
| `emissions_rate` | 11.1 GtC/yr | GCB 2024 |
| Coupling matrix symmetry (MHD module) | — | conserves energy; drift 5.7e-14 |

---

## Current hindcast, 2010–2024

| System | R² | Grade | Note |
|---|---|---|---|
| solar | **+0.359** | C | was −0.443 before the phase fix |
| magnetic | −0.214 | F | finding 18, open |
| atmospheric | +0.471 | C | in-sample |
| oceanic | +0.631 | B | in-sample |
| **mean** | **+0.312** | | |

A mean R² is a weak summary across four incommensurable subsystems — read the
rows, not the average. One subsystem still scores worse than a flat line.

**These are in-sample figures.** For out-of-sample skill see finding 19: held
out on 2020–2024, the shipped parameters score a mean test R² of −0.513.
