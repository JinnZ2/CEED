# Thwaites Glacier — September 2026 risk audit

The [full supplied audit](thwaites_simulators_risk_audit.md) connects Thwaites Glacier studies to cascade, commitment, measurement, and observing-system risk questions. It is preserved verbatim as a contributed research note, including the follow-up on effective redundancy, AMOC coupling, and the THW-F1–THW-F5 falsifier watch list.

## Relationship to CEED

The audit is relevant to [tipping and commitment lag](../../simulation/tipping.py) and [coupled tipping elements](../../simulation/cascade.py). Those links identify existing mechanisms to investigate; this documentation import does not add a calibrated Thwaites model, modify simulation parameters, or wire the audit into the convergence model. Read the [numerical audit](../numerical-audit.md) before interpreting model numbers.

The companion copy is maintained under [`JinnZ2/Simulators/AMOC/studies/thwaites-2026`](https://github.com/JinnZ2/Simulators/tree/docs/thwaites-2026-risk-audit/AMOC/studies/thwaites-2026).

## Supplied materials and provenance

| Material | Role |
|---|---|
| [thwaites_simulators_risk_audit.md](thwaites_simulators_risk_audit.md) | Complete standalone attachment, including the second-pass follow-up dated 2026-09-23. |
| [OKComputer_Thwaites_Glacier_2026_Studies.zip](OKComputer_Thwaites_Glacier_2026_Studies.zip) | Original archive, retained unchanged. It contains one earlier version of the audit, not a collection of paper PDFs. |

The archived Markdown is an exact prefix of the standalone audit. The standalone file appends the follow-up wave; no earlier text was replaced. Both original attachments are retained so that their different contents remain traceable.

| Original attachment | SHA-256 |
|---|---|
| `thwaites_simulators_risk_audit.md` | `080bacf53754558e257aac096e8ebb27b3b9a3df500381b5a3f3276bd3ddb524` |
| `OKComputer_Thwaites_Glacier_2026_Studies.zip` | `1d83974e8ebfbb4de497f2139bfd7be9b62a02441b2ada97aed9e42f5d895555` |

## Evidence boundary

Statements that tools were cloned and run, and the numerical outputs quoted in the audit, belong to the supplied document. This import does not independently reproduce those paper-specific runs or verify the underlying publications. The attachments do not include the custom input files, complete raw run logs, exact simulator commit, or paper PDFs needed for that reproduction.

The chain example uses explicitly synthetic, arbitrary-unit values; it is not a calibrated glacier forecast. The audit also explicitly leaves the Southern-Ocean-to-AMOC transfer path **UNMEASURED**. Its box-model response and declared Sv-to-F calibration must not be treated as measured routing or a real-world collapse threshold.
