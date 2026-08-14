# CEED Numerical Audit

**Scope:** parameters, thresholds, units, and literature-facing values across
the model code and the risk documentation.
**Method:** every internal finding below was computed directly from the
shipped code, not estimated. Literature values were checked against primary
or near-primary sources, cited inline.
**Status:** diagnostic. No model parameters were changed as a result of this
audit — those are scientific choices for the maintainer. Findings are ordered
by severity.

---

## Summary

| # | Finding | Where | Severity |
|---|---------|-------|----------|
| 1 | Phase thresholds are unreachable; every run is Phase 4 from t=0 | `convergence_model.py` | Critical |
| 2 | Retention term makes runaway structural, not emergent | `convergence_model.py` | Critical |
| 3 | Retention-collapse "fix" never engages at realistic temperatures | `minimum_esm_code.py` | Critical |
| 4 | Four incommensurable units are summed into one total | `convergence_model.py` | High |
| 5 | ECS (°C/doubling) used where a forcing (W/m²) is required | `minimum_esm_code.py`, `run_mc.py` | High |
| 6 | Universal model's thresholds sit outside its reachable energy range | `CEED_universal_model.py` | High |
| 7 | Dipole decay rate overstated by ~10x | `red-team-report.md` | High |
| 8 | Roadmap formulas are not dimensionally closed | `To-be-added.md` | Medium |
| 9 | Six declared parameters are never read | `minimum_esm_code.py` | Medium |
| 10 | Emissions rate contradicts its own comment | `minimum_esm_code.py` | Low |
| 11 | Sink split drifts from Global Carbon Budget | `minimum_esm_code.py` | Low |
| 12 | Unsourced scenario probabilities; two overstated observations | `red-team-report.md` | Medium |

Values that check out and should **not** be touched are listed at the end.

---

## 1. Phase thresholds are unreachable — Critical

`ConvergencePredictor` starts at:

```
solar 180.0 + magnetic 92.5 + atmospheric 118.0 + oceanic 110.0 = 500.5
```

The thresholds are:

```
phase_1: 120    phase_2: 150    phase_3: 200    phase_4: 300
```

The initial total exceeds the **Phase 4** threshold by 200.5 units before a
single timestep runs. Over a 3-year run the total energy spans 500.5 → 3021.1,
and `classify_phases` returns:

```
timesteps per phase: {4: 36}
```

All 36 timesteps are Phase 4. Phases 1, 2 and 3 are dead branches — no input
can ever produce them, because energy only grows (finding 2). The phase
classifier currently carries zero information.

`phase_1: 120` is also never compared against anything: `classify_phases`
tests `phase_4`, `phase_3`, `phase_2` and falls through to `1`.

**To make the classifier meaningful**, thresholds need to be on the same scale
as the summed state — roughly 500 (baseline) to 3000 (3-year end state) — or
the state needs normalising to a baseline of 100 per system before summing.

---

## 2. Runaway is structural, not emergent — Critical

Per system the derivative is:

```
dE/dt = input + 1.1*E*(1 - lambda) - lambda*E - 0.001*E^2
```

The retention term `1.1*E*(1-lambda)` is a **level multiplied into a rate**.
Collecting the linear coefficients:

| System | lambda | net linear coeff | implied equilibrium |
|---|---|---|---|
| solar | 0.05 | **+0.9950 /yr** | E ≈ 995 |
| magnetic | 0.02 | **+1.0580 /yr** | E ≈ 1058 |
| atmospheric | 0.08 | **+0.9320 /yr** | E ≈ 932 |
| oceanic | 0.01 | **+1.0790 /yr** | E ≈ 1079 |

Every system has a net linear growth rate near **+100% per year**, checked only
by the `0.001*E²` quadratic. Runaway is therefore not a finding the model
produces — it is imposed by the retention coefficient. The model cannot express
a stable outcome for any input, which makes it unable to answer the question
the project is asking.

The `1.1` is almost certainly meant to be a *retained fraction* (a bit more
energy kept than lost), not a growth multiplier. If retention were written as
`0.1*E*(1-lambda)` — retaining 10% — the linear coefficients would fall to
about +0.045 to +0.089 /yr, and thresholds in the low hundreds would become
meaningful. **This is a modelling decision, not a typo fix, so it is flagged
rather than applied.**

---

## 3. The retention-collapse "fix" never engages — Critical

`minimum_esm_code.py` marks retention collapse as the `CRITICAL FIX` that
prevents unrealistic runaway, and documents it as:

```python
# At T=0: retention ~ 1.05
# At T=3: retention ~ 0.98 (starts losing energy)
# At T=6: retention ~ 0.92 (strong losses)
```

Computed from the shipped `retention_collapse_rate = 0.0008`:

| T (°C) | comment claims | actual | error |
|---|---|---|---|
| 0 | 1.05 | 1.0500 | — |
| 3 | 0.98 | **1.0425** | +0.062 |
| 6 | 0.92 | **1.0202** | +0.100 |

Retention does not cross 1.0 until **T = 7.8 °C**. At every temperature the
model realistically explores, retention stays above 1.0 — that is,
*accumulating*. The stated safeguard is inactive across the entire plausible
range.

To match the documented intent the rate would need to be:

- `0.00767` to hit 0.98 at T=3 (**9.6x** the shipped value)
- `0.00367` to hit 0.92 at T=6 (**4.6x** the shipped value)

Those two targets are mutually inconsistent under a `exp(-k·T²)` form, so the
intended curve cannot be recovered from the comments alone — the maintainer
needs to pick the anchor point. Either the constant is wrong by roughly an
order of magnitude, or the comments are. They cannot both be right.

---

## 4. Incommensurable units are summed — High

`E_current` mixes four different physical quantities and `classify_phases`
adds them:

| Key | Value | Comment says | Actually |
|---|---|---|---|
| `solar` | 180.0 | F10.7 = 180 sfu | solar flux units — genuine |
| `magnetic` | 92.5 | `Kp = 3` | **Kp runs 0–9; 92.5 is not a Kp value** |
| `atmospheric` | 118.0 | — | reads as % of a baseline |
| `oceanic` | 110.0 | — | reads as % of a baseline |

`sum(E)` therefore adds solar flux units to a percentage to a mislabelled
index. The total has no physical meaning, and every threshold comparison
inherits that. Two of the four look like percentages of baseline, which
suggests the intent was a normalised index — in which case `solar: 180` (sfu)
is the outlier and the `Kp = 3` comment is a leftover.

Normalising all four to "% of baseline" would make the sum coherent and would
also put the finding-1 thresholds back in range.

---

## 5. A temperature sensitivity is used as a forcing — High

Both `minimum_esm_code.py` and `run_mc.py` compute:

```python
CO2_forcing = ECS * np.log(CO2 / CO2_preindustrial) / np.log(2)
```

then add the result to `aerosol_ERF` and the solar term, which are in W/m².

But ECS has units of **°C per doubling**, not W/m². At a CO₂ doubling:

- this expression yields **3.00** (which is ECS itself, in °C)
- the standard forcing is `5.35·ln(C/C₀)` = **3.71 W/m²**, and AR6 assesses
  ~3.93 W/m² per doubling

So a temperature is being summed with three forcings and the result fed
through `retention` and `dissipation` as if it were an energy flux. The
numerical coincidence that ECS ≈ 3 and forcing-per-doubling ≈ 3.7–3.9 masks
the error — the magnitudes are close, so the output looks plausible.

The dimensionally correct form separates the two steps: compute forcing
`F = 5.35·ln(C/C₀)` W/m², then convert to temperature via the climate
feedback parameter `λ = F_2x / ECS`.

---

## 6. Universal model thresholds sit outside its reachable range — High

`CEEDSystem` declares `warning=100`, `critical=200`, `tipping=300`, but
`compute_retention` collapses as `1.05·exp(-0.001·E²)`:

| E | retention | dissipation |
|---|---|---|
| 10 | 0.950 | 0.53 |
| 25 | 0.562 | 1.38 |
| 50 | 0.0862 | 2.85 |
| 75 | 0.00379 | 4.40 |
| 100 | **4.77e-05** | 6.00 |
| 300 | **8.60e-40** | 20.20 |

Retention is effectively zero well before the *warning* threshold. Running the
shipped climate example at its shipped forcing:

```
constant forcing     5  -> peak energy    51.89   final state 'stable'
constant forcing    50  -> peak energy    62.86   final state 'stable'
constant forcing   500  -> peak energy   155.11   final state 'stable'
constant forcing  5000  -> peak energy  1551.14   final state 'stable'
```

The demo never exceeds 52 against a warning threshold of 100. Buffer depletion
is gated on `energy > warning_threshold`, so `buffer_capacity` stays pinned at
100% and the `'stressed'`, `'critical'` and `'tipping'` classifications are
unreachable in the shipped example. Roughly 3 orders of magnitude more forcing
is required before any of that machinery activates.

Note also that `compute_retention` and `compute_dissipation` are computed from
`self.state.energy` but the feedback loops are applied to a separate running
`net_effect`, so the two halves of `update()` operate on different quantities.

---

## 7. Dipole decay overstated by ~10x — High

`red-team-report.md` states:

> Earth's dipole field weakening 9–10% per decade in select zones

The measured global rate is about **5% per century**, and the axial dipole is
about **9% weaker than in 1840** — that is, 9% over ~175 years, not per decade
([Olson & Amit](https://pages.jh.edu/polson1/pdfs/ChangesinEarthsDipole.pdf);
[Finlay et al., *Nature Communications*](https://www.nature.com/articles/ncomms10422)).

The "9–10%" figure matches the total decline since 1840 almost exactly, so this
reads as a units slip: a per-175-year change reported as per-decade. As written
the claim is off by roughly a factor of 10 or more.

The South Atlantic Anomaly does weaken faster regionally, which may be the
"select zones" qualifier — but the SAA is also
[not evidence of an imminent reversal](https://www.pnas.org/doi/10.1073/pnas.1722110115),
and the report's framing invites that reading.

---

## 8. Roadmap formulas are not dimensionally closed — Medium

Two formulas in `To-be-added.md` do not resolve to the quantity their name
implies. Both are now implemented verbatim in `simulation/mhd_spatial_model.py`
so the specification is honoured, each paired with a dimensionally correct
companion.

**MHD injection**, `E_input_mhd = η·|v×B|²/ρ`:
`|v×B|²` is (V/m)², dividing by kg/m³ gives V²·m/kg. Not an energy or a flux.
For nominal solar wind (400 km/s, 5 cm⁻³, 5 nT) it returns ~4.8e11 in its own
units. The physical quantities are the Poynting flux (**7.96e-06 W/m²**) and
the bulk kinetic energy flux (**2.68e-04 W/m²**).

**Ionospheric dynamo**, `E_dynamo = ω·(v_atmo×B)·h·σ`:
`|v×B|·σ·h` is a sheet current density (A/m); multiplying by ω (1/s) gives
A/(m·s), not an energy. For the roadmap's own inputs (250 m/s, 45,000 nT,
110 km, 1e-3 S/m) the parts are: sheet current **1.24 A/m**, and Joule
dissipation `σ·|v×B|²·h` = **1.39e-02 W/m²**, which is the physically
meaningful figure and lands in the expected milliwatt band.

Two smaller notes on the roadmap:

- It specifies the dynamo layer at **110 km**. Hall conductivity peaks near
  110 km but Pedersen conductivity — the one that matters for Joule
  dissipation — peaks nearer **125 km**. 110 km is the Hall-weighted choice;
  worth stating explicitly.
- The four spatial buckets are not a partition. Polar / mid-latitude /
  equatorial tile the sphere by latitude (solid-angle fractions 0.134 / 0.366 /
  0.500, summing to 1.000). **Oceanic is a surface type covering ~71% of the
  globe and overlaps all three.** Treating it as a fourth peer band
  double-counts area. It is modelled as a coupled reservoir instead.

The supplied coupling matrix is symmetric with a zero diagonal, which is a good
property: written as `Σ_j C_ij·(E_j − E_i)` it conserves total energy exactly
(measured drift over a 10-year closed run: **5.7e-14**).

---

## 9. Six declared parameters are never read — Medium

`ClimateParameters` declares 17 fields. These 6 are never referenced in
`MinimalESM`:

- `ECS_range` — uncertainty declared, only `ECS_mean` used
- `aerosol_ERF_range` — same
- `GHG_ERF` (3.32) — declared, never applied to any forcing sum
- `land_sink` (0.31), `ocean_sink` (0.23) — only `total_sink_baseline` used
- `permafrost_feedback_range` (14, 175) — `permafrost_feedback()` hardcodes
  `permafrost_mid = 95` instead

The permafrost case is the notable one: 95 is the arithmetic midpoint of
(14, 175), so it is self-consistent, but the ESM ignores the range while
`run_mc.py` correctly samples across it. The two models therefore disagree
about whether permafrost strength is uncertain.

The declared-but-unused ranges make the model look uncertainty-aware in the
deterministic path when it is running fully deterministically.

---

## 10. Emissions rate contradicts its own comment — Low

```python
emissions_rate = 10.0  # GtC/year (current ~11 GtC/year)
```

The comment is closer to right than the value. The Global Carbon Budget 2024
gives total anthropogenic emissions of **11.1 GtC/yr** (2023) and a preliminary
**11.3 GtC/yr** (2024)
([GCB 2024](https://essd.copernicus.org/articles/17/965/2025/)). Using 10.0
understates present-day emissions by ~10%.

At the shipped 54% sink this yields 2.17 ppm/yr growth against an observed
~2.4 ppm/yr, consistent with the understatement.

---

## 11. Sink split drifts from the Global Carbon Budget — Low

| Parameter | Code | GCB 2024 (2014–2023) |
|---|---|---|
| `land_sink` | 0.31 | 0.29 |
| `ocean_sink` | 0.23 | **0.26** |
| `total_sink_baseline` | 0.54 | ~0.55 |

The total is fine; the split is off, with ocean understated by 3 points. Low
impact because only the total is actually used (finding 9) — but if the split
is ever wired in, it should be corrected first.

---

## 12. Risk report: unsourced probabilities and two overstated claims — Medium

**The scenario probability table carries no derivation.** GIC Cascade 40%,
AMOC Stall 25%, Nonlinear Coupling 30%, Atmospheric Arc 10%, Crustal Shear
Slip 5%. These are stated to the nearest 5% with no time horizon, no method,
and no source. A red-team document is strengthened, not weakened, by labelling
these as calibrated analyst judgement with an explicit horizon — as written
they read as model output, and the model (findings 1–4) cannot produce them.

**"Thermosphere still expanded from 2024 storm"** (report dated July 2025).
Storm-time thermospheric density can rise by more than an order of magnitude,
but it recovers over days, not a year
([Parker & Linares 2024](https://arxiv.org/abs/2406.08617)). A residual
expansion persisting 14 months after the Gannon storm is not supported.

**"AMOC slowdown occurring faster than decadal models suggest"**. The RAPID
array shows a weakening trend of **−0.8 ± 0.7 Sv/decade** since 2004 — the
uncertainty is nearly as large as the trend, and the slowdown is not apparent
in the more northerly sections. Tipping-time estimates exist but span wide
intervals (e.g. 2037–2109). "Faster than models suggest" overstates what the
observations currently support; "consistent with the faster end of model
projections, with wide uncertainty" would be defensible
([Carbon Brief AMOC explainer](https://interactive.carbonbrief.org/amoc-explainer/index.html)).

The `phase_2` reference in `analyze_event_timing` compounds finding 1: since
system energy always exceeds `phase_2`, **every** external event is flagged as
critically timed. A 3-year run reports 1,756 of 1,756 events as critical, which
is the same as reporting none.

---

## Values that check out — leave these alone

These were verified and are correct. Do not "fix" them.

| Parameter | Code | Source |
|---|---|---|
| `ECS_mean` = 3.0, likely range (2.5, 4.0) | ✅ | AR6 best estimate 3 °C, likely 2.5–4 °C |
| `aerosol_ERF_mean` = −1.1, range (−1.7, −0.4) | ✅ | AR6 assesses −1.1 W/m², 5–95% −1.7 to −0.4 |
| `permafrost_feedback_range` = (14, 175) GtCO₂/°C | ✅ | [AR6 WG1 FAQ Ch.5](https://www.ipcc.ch/report/ar6/wg1/downloads/faqs/IPCC_AR6_WGI_FAQ_Chapter_05.pdf) |
| `GHG_ERF` = 3.32 W/m² | ✅ | AR6 WMGHG ERF 1750–2019 (declared but unused — finding 9) |
| `total_sink_baseline` = 0.54 | ✅ | GCB ~55% |
| 2.12 GtC = 1 ppm conversion | ✅ | standard |
| `solar_cycle_years` = 11.0 | ✅ | standard |
| F10.7 exceeding 180 sfu | ✅ | plausible for Cycle 25 maximum |
| `CO2_initial` = 420 ppm, `T_initial` = 1.1 °C | ✅ | reasonable for the stated 2024 baseline |
| Coupling matrix symmetry / zero diagonal | ✅ | conserves energy; drift 5.7e-14 |

---

## Suggested order of work

1. **Findings 1–3** decide whether the models can answer their own question.
   Nothing downstream is trustworthy until the retention coefficient and the
   thresholds are on a consistent scale.
2. **Finding 4** (units) should be settled before 1, since normalising the
   state vector likely fixes the thresholds for free.
3. **Finding 5** is a contained change to two functions.
4. **Findings 7, 12** are documentation edits — cheapest to fix, and they are
   the claims a sceptical reader will check first.
5. **Findings 9–11** are hygiene.

---

*Every internal figure in this document was computed from the code as
committed. Rerun the checks with `pytest` and the module `__main__` blocks.*
