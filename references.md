# CEED Scientific References

## Primary Sources

### IPCC AR6 Working Group I (2021)

**“Climate Change 2021: The Physical Science Basis”**

Key parameters used in CEED models:

#### Equilibrium Climate Sensitivity (ECS)

- Best estimate: 3.0°C
- Likely range: 2.5-4.0°C
- Very likely range: 2.0-5.0°C
- Source: AR6 WGI Chapter 7, Executive Summary

#### Aerosol Effective Radiative Forcing

- Total aerosol ERF (1750-2019): -1.1 W/m²
- Very likely range: -1.7 to -0.4 W/m²
- Source: AR6 WGI Chapter 7, Table 7.8

#### Total Greenhouse Gas Forcing

- GHG ERF (excluding ozone/stratospheric H₂O): 3.32 ± 0.29 W/m² (1750-2019)
- Source: AR6 WGI Chapter 7

#### Carbon Sinks

- Land sink: ~31% of anthropogenic CO₂ (2010-2019)
- Ocean sink: ~23% of anthropogenic CO₂ (2010-2019)
- Combined: ~54% total uptake
- Source: AR6 WGI Chapter 5, Executive Summary
- Note: Sink efficiency declines under continued warming

#### Permafrost Carbon Feedback

- Estimated release: 14-175 GtCO₂ per 1°C warming
- High uncertainty due to complex thaw dynamics
- Source: AR6 WGI Chapter 5, FAQ 5.1
- Mid-range estimate used in model: ~95 GtCO₂ per °C

-----

## Additional Literature

### Solar Cycle Variability

- 11-year solar cycle
- Total Solar Irradiance (TSI) variation: ~0.1%
- Radiative forcing variation: ~0.1 W/m²
- Source: Kopp & Lean (2011), “A new, lower value of total solar irradiance”

### Cloud Feedback

- Remains largest source of uncertainty in climate sensitivity
- Net positive feedback in most models
- Source: Sherwood et al. (2020), “An Assessment of Earth’s Climate Sensitivity”

### Retention and Dissipation Mechanisms

- High-energy plasma states exhibit increased leakage
- Radiative cooling scales approximately as T⁴ (Stefan-Boltzmann)
- Simplified as T^1.2 in minimal model for computational efficiency

-----

## Model Simplifications and Assumptions

### What We Include

✓ Key feedback mechanisms (water vapor, ice-albedo, radiative cooling, carbon sinks)
✓ Solar cycle integration
✓ Aerosol forcing scenarios
✓ Permafrost carbon feedback
✓ Retention collapse at high energy (prevents unrealistic runaway)

### What We Simplify

⚠️ Single-layer atmosphere (real climate has complex vertical structure)
⚠️ Lumped ocean (no deep ocean representation)
⚠️ Simplified carbon cycle (no explicit biosphere dynamics)
⚠️ Idealized feedback forms (real feedbacks have complex spatial patterns)

### What We Omit

❌ Spatial dimensions (this is a 0-D energy balance model)
❌ Ice sheet dynamics (multi-century timescales)
❌ Detailed ocean circulation
❌ Land use change
❌ Non-CO₂ greenhouse gases (implicitly included in GHG forcing)

-----

## Validation Approach

Our minimal model is designed to:

1. Capture essential feedback dynamics
1. Use IPCC-constrained parameter ranges
1. Produce realistic temperature trajectories
1. Enable rapid scenario exploration

**It is NOT:**

- A substitute for comprehensive Earth System Models
- Suitable for regional or spatial predictions
- Capable of capturing all physical processes

**It IS useful for:**

- Understanding feedback interactions
- Exploring parameter sensitivity
- Testing policy scenarios
- Educational purposes
- Rapid prototyping of ideas

-----

## Comparison to Other Models

### Complexity Spectrum

**Simple Energy Balance Models (EBMs):**

- 0-1 spatial dimensions
- Fast computation
- Limited process representation
- Our model fits here

**Intermediate Complexity Models:**

- 2-3 spatial dimensions
- Some process detail
- Hours to days for multi-century runs

**Comprehensive ESMs (Earth System Models):**

- 3D atmosphere, ocean, land, ice
- Detailed process representation
- Weeks to months for multi-century runs
- IPCC projections use these

-----

## Key Papers for Deeper Reading

### Feedback Dynamics

- Roe & Baker (2007): “Why Is Climate Sensitivity So Unpredictable?”
- Knutti et al. (2017): “Beyond equilibrium climate sensitivity”

### Tipping Points

- Lenton et al. (2008): “Tipping elements in the Earth’s climate system”
- Armstrong McKay et al. (2022): “Exceeding 1.5°C global warming could trigger multiple climate tipping points”

### Carbon Cycle

- Friedlingstein et al. (2020): “Global Carbon Budget 2020”
- IPCC AR6 WGI Chapter 5: “Global Carbon and other Biogeochemical Cycles and Feedbacks”

### Aerosols

- Bellouin et al. (2020): “Bounding global aerosol radiative forcing of climate change”
- IPCC AR6 WGI Chapter 7: “The Earth’s Energy Budget, Climate Feedbacks, and Climate Sensitivity”

-----

## Data Sources for Real-Time Integration

### When moving from mock data to real inputs:

**Solar Activity:**

- NOAA Space Weather Prediction Center: https://www.swpc.noaa.gov/
- F10.7 solar flux index (daily updates)

**Geomagnetic Activity:**

- Kp index: https://www.gfz-potsdam.de/en/kp-index/
- NOAA Space Weather: https://www.swpc.noaa.gov/products/planetary-k-index

**Atmospheric Data:**

- NOAA Global Monitoring Laboratory: https://gml.noaa.gov/
- CO₂ concentration at Mauna Loa

**Ocean Circulation:**

- RAPID-AMOC Array: https://rapid.ac.uk/
- Atlantic Meridional Overturning Circulation monitoring

-----

## Updates and Corrections

This is a living document. As new research emerges or model refinements occur, references will be updated.

**Last updated:** [Date]
**CEED Version:** 1.0

-----

## Citation

If you use or build upon CEED, please cite:

```
CEED: Cascading Energy Event Disruption Framework
[Repository URL]
A minimal Earth system model and universal feedback framework
Co-created by [contributors]
2025
```

-----

## Contact for Scientific Questions

Open an issue on the GitHub repository for:

- Parameter questions
- Literature suggestions
- Model improvements
- Scientific discussion

We welcome contributions from domain experts!
