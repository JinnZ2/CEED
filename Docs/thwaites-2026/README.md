# Thwaites Glacier — September 2026 risk audit

The [full supplied audit](thwaites_simulators_risk_audit.md) connects Thwaites Glacier studies to cascade, commitment, measurement, and observing-system risk questions. It is preserved verbatim as a contributed research note, including the follow-up on effective redundancy, AMOC coupling, the THW-F1–THW-F5 falsifier watch list, and section 9 on ENSO forcing pulses.

## Relationship to CEED

The audit is relevant to [tipping and commitment lag](../../simulation/tipping.py) and [coupled tipping elements](../../simulation/cascade.py). Those links identify existing mechanisms to investigate; this documentation import does not add a calibrated Thwaites model, modify simulation parameters, or wire the audit into the convergence model. Read the [numerical audit](../numerical-audit.md) before interpreting model numbers.

The companion copy is maintained under [`JinnZ2/Simulators/AMOC/studies/thwaites-2026`](https://github.com/JinnZ2/Simulators/tree/main/AMOC/studies/thwaites-2026).

## Supplied materials and provenance

| Material | Role |
|---|---|
| [thwaites_simulators_risk_audit.md](thwaites_simulators_risk_audit.md) | Revised standalone attachment, including the second-pass follow-up and added section 9 on ENSO coupling. |
| [OKComputer_Thwaites_Glacier_2026_Studies.zip](OKComputer_Thwaites_Glacier_2026_Studies.zip) | Original archive, retained unchanged. It contains one earlier version of the audit, not a collection of paper PDFs. |

The current audit reproduces the revised attachment `thwaites_simulators_risk_audit(1).md` byte-for-byte under the canonical filename `thwaites_simulators_risk_audit.md`. Both the archived Markdown and the previous standalone attachment are exact prefixes of this revision: section 9 is appended without changing earlier text. The ZIP is unchanged, and the previous standalone version remains in Git history.

| Supplied version | SHA-256 |
|---|---|
| Current audit, supplied as `thwaites_simulators_risk_audit(1).md` | `215822455536a97ddc277d886a0f5f0fd3668e6f3014fa11124b35f9664ab124` |
| Previous standalone audit (Git history) | `080bacf53754558e257aac096e8ebb27b3b9a3df500381b5a3f3276bd3ddb524` |
| `OKComputer_Thwaites_Glacier_2026_Studies.zip` | `1d83974e8ebfbb4de497f2139bfd7be9b62a02441b2ada97aed9e42f5d895555` |

## Evidence boundary

Statements that tools were cloned and run, and the numerical outputs quoted in the audit, belong to the supplied document. This import does not independently reproduce those paper-specific runs or verify the underlying publications. The attachments do not include the custom input files, complete raw run logs, exact simulator commit, or paper PDFs needed for that reproduction.

The chain example uses explicitly synthetic, arbitrary-unit values; it is not a calibrated glacier forecast. The audit also explicitly leaves the Southern-Ocean-to-AMOC transfer path **UNMEASURED**. Its box-model response and declared Sv-to-F calibration must not be treated as measured routing or a real-world collapse threshold.

Section 9 carries the source date **2026-09-24**, preserved as supplied. Its time-sensitive ENSO statements and illustrative melt arithmetic have not been independently validated in this import; the source itself labels the arithmetic as declared, not a result.

## Further reading to investigate

The [multilingual citation backlog](CITATION_BACKLOG.md) adds three September 2026 papers and three older 2026 references with specific follow-up questions. All six papers are English-language; non-English discovery sources are labeled separately. No qualifying original non-English paper was verified in the bounded 1–23 September search. These citations are a reading queue, not completed simulator analyses.

## Note 2026-10-28: a second freshwater term since 2016 (co-authored)

Mishra, Gayen & Naveira Garabato 2026, *Reversal of Antarctic Intermediate Water trends triggered by sea ice decline*, Nature Climate Change, doi:10.1038/s41558-026-02758-y. Read from press coverage only; the paper was not reachable.

| | Freshwater term | Direction since 2016 | Water mass it feeds |
|---|---|---|---|
| This audit, CDW-AABW loop | ice-shelf basal melt (Thwaites, Amundsen Sea) | up | AABW, continental shelf |
| AAIW study | summer sea-ice meltwater, about 3x -> 2x the Amazon's annual flow | down, about 6,600 km³/yr (0.21 Sv) | AAIW, toward the polar front |

The sea-ice change is about 63x the CDW-attributable basal-melt term in `JinnZ2/earth-systems-physics` `layer_4_hydrosphere.py` (104.8 Gt/yr), computed in that repo's `aaiw_reversal_2026.py` (DERIVED). PROPOSED reading: a Southern Ocean freshwater budget that carries only the ice-shelf term is missing the larger moving term, if the stated melt volumes are net delivery rather than gross seasonal melt (VERIFY in the paper). The coverage does not mention ice-shelf meltwater, so whether the AAIW budget carries Thwaites is OWED.

Both studies date their change from the 2015-16 sea-ice drop; Spira et al. 2026 (Nature Climate Change, doi:10.1038/s41558-026-02601-4) attribute that drop to wind-driven mixing over a thinned Winter Water layer. Nothing here is wired into CEED's models.
