# Ionospheric Chaos Index (ICI) Module
*Extension to model thermospheric and ionospheric instability from multi-system and lunar perigee forcing.*

---

##  Purpose

Capture the chaotic, non-periodic fluctuations in Earth's upper atmosphere due to:
- Lunar perigee and declination
- Solar wind enhancements
- Geomagnetic field weakening
- Ring current persistence
- Gravity wave amplification

---

##  Core Idea

**Ionosphere is not just a passive recipient** of solar radiation—it's a **coupled oscillator** influenced by:
- **Lunar tidal harmonics**
- **Thermal instabilities**
- **Magnetospheric charge coupling**

Use this module to create an **instability metric**.

---

##  ICI Equation

```python
ICI(t) = α_solar(t) · ρ_density(t) · Λ_geomag(t) · ξ_lunar(t) + ε_noise(t)

Where:
	•	α_solar(t) = solar forcing index (e.g. F10.7 flux)
	•	ρ_density(t) = thermospheric density variation (from CEED)
	•	Λ_geomag(t) = magnetospheric decay factor
	•	ξ_lunar(t) = normalized lunar perigee contribution
	•	ε_noise(t) = stochastic turbulence term (high-pass filtered)

⸻

 CEED Integration
	•	Use ICI(t) as a modulator of atmospheric energy retention
	•	Inject into resonance_coupling() to simulate chaotic phase transfer
	•	Add to stochastic_forcing() to reflect upper-atmosphere disorder

⸻

 Lunar Coupling Input (Optional)

ξ_lunar(t) = 1 + 0.1 * sin(2π * t / 27.3) + 0.05 * sin(2π * t / 206)

Where:
	•	27.3 days = sidereal month (orbital)
	•	206 days = declination return (nodding)

Impact

High ICI correlates with:
	•	Radio blackouts
	•	Aurora anomalies
	•	Satellite position drift
	•	Amplified thermal upwelling

When CEED predicts ICI spikes, it’s time to worry about tech failures and cross-system chaos triggers.

⸻

 Supporting Concepts
	•	Lunar Ionospheric Forcing (LIF)
	•	Heliospheric current sheet instability
	•	D-region conductivity inversions
	•	Auroral heating feedback loops

⸻

 Suggested Outputs
	•	ICI plot over time
	•	ICI → probability of GPS drift event
	•	ICI thresholds to trigger tech warning layer in CEED dashboard

