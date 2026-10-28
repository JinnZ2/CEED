# Arase, May 2024 superstorm: multi-MeV electrons inside the "impenetrable barrier"

Note 2026-10-28: co-authored, human + Claude. Documentation import only. No simulation code, parameter or threshold is changed.

## The study

Miyoshi et al. 2026, *Local Acceleration of Multi-MeV Electrons Inside Earth's Impenetrable Barrier During the May 2024 Superstorm*, Geophysical Research Letters, doi:10.1029/2026GL125879.

| Read | Not read |
|---|---|
| [Nagoya University press release](https://en.nagoya-u.ac.jp/news/articles/pr-2024-solar-storm-triggered-ultra-relativistic-electrons-inside-earths-impenetrable-barrier/), 2026-10-28 | the paper, figures, methods, supplement (Wiley and phys.org returned client errors through the fetch tool) |
| [CALET on the ISS, arXiv:2602.03990](https://arxiv.org/abs/2602.03990), abstract page: independent instrument, same storm | Baker et al. 2014, where the barrier was first reported |

Full record: [`JinnZ2/study-aperture/examples/arase-may2024`](https://github.com/JinnZ2/study-aperture/tree/arase-may2024/examples/arase-may2024) (PR #8, open).

## What was found, at press-release depth

| Claim | Status here |
|---|---|
| Electrons above 3 MeV were present inside about 3 Earth radii after the storm | found: Arase and CALET independently (CALET: L 2.2-3.2, below the nominal L 2.8) |
| They were accelerated locally, not carried in from the outer belt | unresolved: rests on a growing phase-space-density peak whose coordinates depend on a magnetic field model not named in the material read; loss is not carried in the inference |
| Chorus waves did it, helped by an eroded plasmasphere (recovery about 4.3-4.7 days vs about 2 typical) | unresolved: timing only; the release says "potentially" |
| The barrier was "broken" | the barrier is an empirical boundary from observations up to 2014; a superstorm is outside the conditions it was characterized under (PROPOSED reading) |

CALET lifetimes as read (VERIFY against the paper): >1.5 MeV persisted over 5 months; e-folding about 41 / 25 / 19 days for >1.5 / >3.4 / >8.2 MeV at L 2.8-3.2; about 18 / 45 / 72 days at L 2.3.

## Relationship to CEED

| CEED piece | Bearing |
|---|---|
| [red-team report](../../red-team-report.md), "Van Allen belts not dissipating normally" | CALET measures months-long persistence of a new low-L population after May 2024. Whether that is "not normal" needs a stated baseline of decay times at the same L and energy; none is given in either document. See the dated note beside that line. |
| [`simulation/tipping.py`](../../simulation/tipping.py), event-induced tipping | Same shape (PROPOSED): a boundary that holds under the typical range is crossed by one large excursion. A boundary characterized on moderate storms is a reading of those storms, not a wall. |
| Magnetic subsystem (Kp proxy, held-out test R² +0.057, finding 18 in [numerical-audit](../numerical-audit.md)) | Kp is a global index; it cannot carry where in L energy is deposited. This event is a case the current magnetic representation has no slot for. Out of scope to wire in here; maintainer decision. |

## Evidence boundary

Nothing here verifies the paper's method. The mechanism question stays open until the phase-space-density peak is checked under a second field model and against a transport-plus-loss alternative.
