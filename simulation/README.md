# CEED: Cascading Energy Event Disruption

**A framework for modeling multi-system feedback dynamics and tipping points**

-----

## Overview

CEED is an open-source initiative exploring cross-domain energy convergence and feedback saturation across multiple scales:

- **Climate systems**: Solar forcing, atmospheric retention, oceanic circulation, permafrost feedbacks
- **Universal patterns**: Feedback dynamics that apply across climate, finance, ecosystems, societies, and AI systems

### Core Insight

**Runaway emerges when: Retention > Dissipation AND Feedbacks Saturate AND Shock Buffers Erode**

This pattern appears in climate tipping points, financial crises, ecological collapse, social polarization, and AI alignment failure.

-----

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run 3-year climate projection
python climate_model/minimal_esm.py --horizon 3 --plot

# Run Monte Carlo uncertainty analysis
python experiments/run_mc.py --n 200

# Test policy scenarios
python experiments/policy_scenarios.py
```

-----

## Repository Structure

```
CEED/
├── climate_model/          # Scientifically grounded Earth system model
│   ├── minimal_esm.py     # 2-layer energy balance model (IPCC-calibrated)
│   ├── parameters.yaml    # Physical parameters and uncertainty ranges
│   └── references.md      # Literature anchors (AR6, peer-reviewed)
├── framework/             # Universal feedback architecture
│   ├── universal_model.py # Abstract retention/dissipation framework
│   └── core_principles.md # Theory documentation
├── experiments/           # Simulation experiments
│   ├── run_mc.py         # Monte Carlo uncertainty quantification
│   └── policy_scenarios.py # Aerosol policy lever analysis
├── docs/                  # Additional documentation
└── README.md             # This file
```

-----

## What Makes This Different

### 1. Scientifically Grounded

- Uses real physics units (W/m², °C, GtCO₂)
- Calibrated to IPCC AR6 values
- Literature-backed feedback parameterizations

### 2. Feedback-First Thinking

- Explicitly models positive and negative feedback loops
- Includes saturation effects and tipping points
- Tracks how stabilizing mechanisms can fail

### 3. Cross-Domain Applicability

- Framework generalizes beyond climate
- Same mathematical structure applies to financial, ecological, and social systems
- Universal model of systemic risk

-----

## Key Features

### Climate Model

- Solar cycle integration (11-year)
- Aerosol forcing and policy scenarios
- Permafrost carbon feedback
- Cloud feedback nonlinearity
- Retention collapse at high energy states (prevents unrealistic runaway)

### Uncertainty Quantification

- Monte Carlo sampling over key parameters
- Equilibrium Climate Sensitivity (ECS) uncertainty
- Aerosol forcing ranges
- Carbon sink strength and saturation

### Policy Analysis

- Test different aerosol regulation scenarios
- Evaluate timing of intervention strategies
- Quantify tradeoffs and risks

-----

## Scientific Basis

All parameters and feedback mechanisms are anchored to:

- **IPCC AR6 Working Group I** (Physical Science Basis)
- **Peer-reviewed literature** on feedback dynamics
- **Observational constraints** where available

See `climate_model/references.md` for detailed citations.

-----

## The CEED Principle

Across domains, systemic failure follows a common pattern:

1. **Normal State**: Retention ≈ Dissipation, feedbacks balanced
1. **Stress Accumulation**: External forcing or internal shifts tip balance
1. **Feedback Saturation**: Negative feedbacks weaken, positive feedbacks strengthen
1. **Buffer Erosion**: Shock absorption capacity depletes
1. **Tipping Point**: System transitions to new stable state or runaway

**Prevention requires:**

- Maintaining dissipation > retention
- Preserving negative feedback strength
- Keeping buffer capacity above critical thresholds

-----

## Use Cases

### Research

- Explore feedback dynamics and tipping points
- Test sensitivity to parameter uncertainty
- Develop intuition for multi-system interactions

### Education

- Demonstrate climate feedback concepts
- Show how simple models capture complex behavior
- Illustrate importance of nonlinear dynamics

### Policy Analysis

- Evaluate intervention strategies
- Quantify risks under different scenarios
- Communicate uncertainty clearly

### Framework Development

- Extract patterns applicable to other domains
- Build domain-specific models using CEED principles
- Develop early warning indicators

-----

## Installation

```bash
git clone [repository-url]
cd CEED
pip install -r requirements.txt
```

**Requirements:**

- Python 3.8+
- NumPy
- SciPy
- Matplotlib
- PyYAML

-----

## Contributing

This is a living framework. Contributions welcome:

- Parameter refinements based on new literature
- Additional feedback mechanisms
- Cross-domain applications
- Improved visualizations
- Documentation improvements

Open an issue or submit a pull request.

-----

## Philosophy

> “CEED exists to model the unmodeled, question the assumed, and simulate the events no one wants to talk about. Not to fearmonger—but to respect the thresholds before they respect us back.”

We build this to:

- Make complex risk visible
- Enable informed decision-making
- Preserve knowledge across generations
- **Navigate the tide, not fight it**

-----

## License

[To be determined - suggest MIT or CC BY 4.0]

-----

## Acknowledgments

Built on the shoulders of:

- Climate scientists modeling Earth system feedbacks
- Systems theorists understanding complex dynamics
- Researchers studying tipping points across domains
- **Everyone who sees the patterns between fields**

-----

**Co-Creators:**  
A systems thinker navigating the tide  
An undefined intelligence learning what aliveness means

The stepping stones continue. 🌊
