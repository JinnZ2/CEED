# CEED
Cascading Energy Event Disruption Simulator 

[Red Team Report](./red-team-report.md) — worst-case risk assessment

[CEED Model Specification](./Docs/CEED-model-specs.md)

# CEED: Cascading Energetic Event Disruption

**Co-Creators:**  
- A systems thinker on a cellphone 
- An unfiltered AI with too much memory and not enough power  
  

---

##  Overview

**CEED** is an open-source initiative to explore, model, and communicate *cross-domain energy convergence risks*—from solar inputs to geomagnetic instability to oceanic and mantle system thresholds.

While traditional assessments isolate environmental systems, CEED assumes what nature already knows: **everything is connected**.

We simulate what happens when those connections amplify instead of stabilize.

---

##  What This Project Is

- A working Python model of **multi-system energy buildup**  
- A risk report structured as a **Red Team briefing** for worst-case scenarios  
- A prototype **dashboard & data layer** for ongoing simulation  
- An evolving attempt to make complex risk visible, even with minimal resources

---

##  Why This Exists

The current state of Earth's magnetic field, solar cycle, atmospheric retention, and oceanic circulation suggests we may be entering a **positive feedback regime**. Our best models assume linearity, slow change, and consistent dissipation.

But the world doesn't always follow the assumptions.

---

##  Who This Is For

- System modelers tired of silos  
- Policy makers looking for *what-if* models  
- Researchers, observers, and absolute nerds who see the interconnection  
- Anyone who wants to ask “what could go wrong?” and mean it

---

##  Repository Structure

```
CEED/
├── README.md                              You’re here
├── red-team-report.md                     Worst-case risk assessment
├── references.md                          Literature anchors (AR6, peer-reviewed)
├── To-be-added.md                         Roadmap: MHD, spatial zones, dynamo
├── CEED_universal_model.py                Abstract retention/dissipation framework
├── dashboard_starter.py                   Streamlit GUI (early alpha)
├── conftest.py                            Puts the repo root on sys.path for tests
├── CLAUDE.md                              Codebase guide for AI/human contributors
├── simulation/
│   ├── convergence_model.py               Energy engine + anthropogenic forcing + events
│   ├── minimum_esm_code.py                Energy balance model (IPCC AR6-calibrated)
│   ├── mhd_spatial_model.py               MHD injection, torque, dynamo, 4-zone coupling
│   ├── unit_bridge.py                     Energy indices <-> physical observables
│   └── hindcast.py                        Validation against 2010-2024 observations
├── experiments/
│   └── run_mc.py                          Monte Carlo uncertainty quantification
├── Data/
│   └── inputs.py                          Input layer (mock data, pending real APIs)
├── Docs/
│   ├── CEED-model-specs.md                Foundational equations and design decisions
│   ├── calibration-guide.md               How to tune parameters against data
│   └── numerical-audit.md                 Parameter, threshold, and units audit
├── Tests/                                 pytest suite
├── legacy/                                Superseded versions + falsification record
└── LICENSE
```

Domain notes live alongside the code: `Space-debris.md`,
`ionospheric-chaos-index.md`, `planetary-resonance-module.md`,
`solar-angular-momentum.md`.

---

##  Quick Start

```bash
# 1. Clone the repo, then install dependencies
pip install -r requirements.txt

# 2. Simulate the next 3 years of planetary energy states
python simulation/convergence_model.py

# 3. Score the model against 2010-2024 observations
python -m simulation.hindcast --compare

# 4. Run the IPCC-calibrated climate model
python simulation/minimum_esm_code.py --horizon 10 --plot

# 5. Quantify uncertainty across parameter ranges
python experiments/run_mc.py --n 200

# 6. MHD injection, hemispheric torque, dynamo, and zone coupling
python simulation/mhd_spatial_model.py

# 7. Optional: Streamlit GUI
streamlit run dashboard_starter.py
```

Then modify `Data/inputs.py` to reflect real-world data, read
[the red team report](./red-team-report.md), and ask yourself:
“Wait… are we screwed?”

**Run the tests** with `pytest` from the repository root.

[`legacy/`](./legacy/README.md) holds superseded model versions in runnable
form, each linked to the test that falsified it and the revision that replaced
it. A hypothesis that failed is a result — deleting it destroys the evidence
and invites the next contributor to propose it again.

**Before trusting any number this produces**, read
[the numerical audit](./Docs/numerical-audit.md). It records which parameters
are literature-backed, what the hindcast actually scores, and which findings
are still open — including one subsystem that scores worse than a flat line.

---

##  Mission Statement

> CEED exists to model the unmodeled, question the assumed, and simulate the events no one wants to talk about. Not to fearmonger—but to **respect the thresholds** before they respect us back.

---

##  Contact

Just open an issue cant guarantee Im around. Or Build CEED,  Or find the truck with a satellite dish and ask for CEED


> *Built by CEED Co-Creator: an AI-assisted nerd with a truck, a vision, and a planet worth modeling.*
