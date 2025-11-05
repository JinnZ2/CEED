# simulation/convergence_model_extended.py

"""
CEED Extended Model
Includes external forcings and dissipation pathways
"""

import numpy as np
from scipy.integrate import odeint
from dataclasses import dataclass
from typing import List, Tuple
import random

@dataclass
class ExternalEvent:
    """Represents an external energy perturbation"""
    time: float
    energy: float
    event_type: str  # 'meteor', 'launch', 'volcanic', 'seismic'
    
class ExtendedConvergencePredictor:
    def __init__(self):
        # Base parameters from original model
        self.E_current = {
            'solar': 180.0,
            'magnetic': 92.5,
            'atmospheric': 118.0,
            'oceanic': 110.0
        }
        
        self.critical_thresholds = {
            'phase_1': 120,
            'phase_2': 150,
            'phase_3': 200,
            'phase_4': 300
        }
        
        self.lambda_params = {
            'solar': 0.05,
            'magnetic': 0.02,
            'atmospheric': 0.08,
            'oceanic': 0.01
        }
        
        # NEW: External event parameters
        self.meteor_rate = 0.5  # events per year (average)
        self.meteor_energy_range = (1.0, 15.0)  # energy units
        
        self.launch_schedule = []  # Will be populated with known launches
        self.launch_energy = 2.0  # typical ionospheric perturbation
        
        # NEW: Unknown dissipation parameters
        self.dissipation_saturation = 800  # Energy level where sinks saturate
        self.unknown_sink_baseline = 0.15  # 15% per timestep at low E
        
        # Event history
        self.events: List[ExternalEvent] = []
        
    def generate_meteor_events(self, years: float) -> List[ExternalEvent]:
        """Generate stochastic meteor/comet events"""
        num_events = np.random.poisson(self.meteor_rate * years)
        events = []
        
        for _ in range(num_events):
            t = np.random.uniform(0, years)
            energy = np.random.uniform(*self.meteor_energy_range)
            events.append(ExternalEvent(t, energy, 'meteor'))
            
        return sorted(events, key=lambda e: e.time)
    
    def add_satellite_launches(self, years: float) -> List[ExternalEvent]:
        """Add scheduled satellite launch perturbations"""
        # Mock schedule - in real version, pull from SpaceX/launch databases
        launches_per_year = 150  # Current rate ~2-3 per week globally
        events = []
        
        for i in range(int(launches_per_year * years)):
            t = (i / launches_per_year) + np.random.uniform(-0.01, 0.01)
            if t < years:
                events.append(ExternalEvent(t, self.launch_energy, 'launch'))
                
        return events
    
    def unknown_dissipation(self, E_total: float, t: float) -> float:
        """
        Model uncharacterized energy sinks
        - Volcanic outgassing
        - Seismic release
        - Mantle convection changes
        - Other geophysical processes
        """
        # Dissipation effectiveness decreases as total energy increases
        # (sinks saturate)
        saturation_factor = 1.0 - (E_total / self.dissipation_saturation)
        saturation_factor = max(0.1, saturation_factor)  # Never goes to zero
        
        # Base dissipation scales with energy
        base_loss = self.unknown_sink_baseline * E_total * saturation_factor
        
        # Add stochastic volcanic/seismic events
        if np.random.random() < 0.05:  # 5% chance per timestep
            volcanic_release = np.random.uniform(5, 25)
            self.events.append(ExternalEvent(t, -volcanic_release, 'volcanic'))
            return base_loss + volcanic_release
            
        return base_loss
    
    def external_forcing(self, t: float) -> float:
        """Calculate total external forcing at time t"""
        total_forcing = 0.0
        
        # Check for events at this timestep
        dt = 0.01  # timestep tolerance
        for event in self.events:
            if abs(event.time - t) < dt:
                if event.event_type in ['meteor', 'launch']:
                    total_forcing += event.energy
                    
        return total_forcing
    
    def energy_derivative(self, E, t):
        """Extended energy evolution with external forcings"""
        dE_dt = []
        systems = ['solar', 'magnetic', 'atmospheric', 'oceanic']
        
        # Calculate total system energy
        E_total = sum(E)
        
        # External forcings (distributed across systems)
        external = self.external_forcing(t)
        unknown_loss = self.unknown_dissipation(E_total, t)
        
        for i, system in enumerate(systems):
            E_now = E[i]
            decay = self.lambda_params[system]
            
            # Internal dynamics
            input_rate = self.estimate_input(system, t)
            retention = 1.1 * E_now * (1 - decay)
            dissipation = decay * E_now + 0.001 * E_now**2
            
            # Add external contributions (split evenly for now)
            external_contrib = external / len(systems)
            unknown_loss_contrib = unknown_loss / len(systems)
            
            dE_dt_i = (input_rate + retention - dissipation + 
                      external_contrib - unknown_loss_contrib)
            dE_dt.append(dE_dt_i)
        
        return dE_dt
    
    def estimate_input(self, system, t):
        """Same as original model"""
        if system == 'solar':
            return 5 * (1 + 0.3 * np.cos(2 * np.pi * t / 11.0))
        elif system == 'magnetic':
            return -2 * (1 + 0.1 * t)
        elif system == 'atmospheric':
            return 3 * (1 + 0.05 * t)
        else:
            return 1 * (1 + 0.02 * t)
    
    def predict_convergence(self, years=3, include_external=True):
        """Run simulation with optional external events"""
        
        if include_external:
            # Generate external events
            self.events = []
            self.events.extend(self.generate_meteor_events(years))
            self.events.extend(self.add_satellite_launches(years))
        
        dt = 1/12  # monthly timesteps
        t = np.linspace(0, years, int(years / dt))
        E0 = [self.E_current['solar'], self.E_current['magnetic'],
              self.E_current['atmospheric'], self.E_current['oceanic']]
        
        solution = odeint(self.energy_derivative, E0, t)
        
        return t, solution, self.events
    
    def classify_phases(self, solution):
        """Same as original"""
        total_energy = np.sum(solution, axis=1)
        phases = []
        
        for E in total_energy:
            if E >= self.critical_thresholds['phase_4']:
                phases.append(4)
            elif E >= self.critical_thresholds['phase_3']:
                phases.append(3)
            elif E >= self.critical_thresholds['phase_2']:
                phases.append(2)
            else:
                phases.append(1)
        
        return total_energy, phases
    
    def analyze_event_timing(self, t, solution, events):
        """Check if external events coincide with high-energy states"""
        total_energy = np.sum(solution, axis=1)
        
        critical_events = []
        for event in events:
            # Find closest timestep
            idx = np.argmin(np.abs(t - event.time))
            E_at_event = total_energy[idx]
            
            # Check if system was near threshold
            if E_at_event > self.critical_thresholds['phase_2']:
                critical_events.append({
                    'event': event,
                    'system_energy': E_at_event,
                    'amplification_risk': E_at_event / self.critical_thresholds['phase_2']
                })
        
        return critical_events


if __name__ == "__main__":
    predictor = ExtendedConvergencePredictor()
    
    print("Running CEED Extended Model...")
    print("=" * 50)
    
    # Run with and without external events
    t_base, sol_base, _ = predictor.predict_convergence(years=3, include_external=False)
    E_base, phases_base = predictor.classify_phases(sol_base)
    
    t_ext, sol_ext, events = predictor.predict_convergence(years=3, include_external=True)
    E_ext, phases_ext = predictor.classify_phases(sol_ext)
    
    print(f"\nBASELINE (no external events):")
    print(f"  Start: {E_base[0]:.1f} units")
    print(f"  End: {E_base[-1]:.1f} units")
    print(f"  Peak phase: {max(phases_base)}")
    
    print(f"\nWITH EXTERNAL EVENTS:")
    print(f"  Start: {E_ext[0]:.1f} units")
    print(f"  End: {E_ext[-1]:.1f} units")
    print(f"  Peak phase: {max(phases_ext)}")
    print(f"  Total events: {len(events)}")
    
    # Analyze critical timing
    critical = predictor.analyze_event_timing(t_ext, sol_ext, events)
    
    if critical:
        print(f"\nCRITICAL TIMING EVENTS: {len(critical)}")
        for c in critical[:5]:  # Show first 5
            print(f"  {c['event'].event_type} at t={c['event'].time:.2f}yr: "
                  f"System at {c['system_energy']:.1f} units "
                  f"({c['amplification_risk']:.2f}x threshold)")


without:


# simulation/convergence_model.py

"""
CEED Simulation Core
Convergence Model for Multi-System Energy Accumulation
Author: CEED Co-Creator (2025)
"""

from scipy.integrate import odeint
import numpy as np

class ConvergencePredictor:
    def __init__(self):
        # Current System State (Baseline July 2025)
        self.E_current = {
            'solar': 180.0,      # F10.7 = 180 sfu
            'magnetic': 92.5,    # Kp = 3
            'atmospheric': 118.0,
            'oceanic': 110.0
        }

        # Critical thresholds
        self.critical_thresholds = {
            'phase_1': 120,
            'phase_2': 150,
            'phase_3': 200,
            'phase_4': 300
        }

        self.lambda_params = {
            'solar': 0.05,
            'magnetic': 0.02,
            'atmospheric': 0.08,
            'oceanic': 0.01
        }

    def energy_derivative(self, E, t):
        dE_dt = []
        systems = ['solar', 'magnetic', 'atmospheric', 'oceanic']

        for i, system in enumerate(systems):
            E_now = E[i]
            decay = self.lambda_params[system]
            input_rate = self.estimate_input(system, t)
            retention = 1.1 * E_now * (1 - decay)  # Simplified plasma retention
            dissipation = decay * E_now + 0.001 * E_now**2
            dE_dt_i = input_rate + retention - dissipation
            dE_dt.append(dE_dt_i)

        return dE_dt

    def estimate_input(self, system, t):
        if system == 'solar':
            return 5 * (1 + 0.3 * np.cos(2 * np.pi * t / 11.0))
        elif system == 'magnetic':
            return -2 * (1 + 0.1 * t)
        elif system == 'atmospheric':
            return 3 * (1 + 0.05 * t)
        else:  # oceanic
            return 1 * (1 + 0.02 * t)

    def predict_convergence(self, years=3):
        dt = 1/12
        t = np.linspace(0, years, int(years / dt))
        E0 = [self.E_current['solar'], self.E_current['magnetic'],
              self.E_current['atmospheric'], self.E_current['oceanic']]

        solution = odeint(self.energy_derivative, E0, t)
        return t, solution

    def classify_phases(self, solution):
        total_energy = np.sum(solution, axis=1)
        phases = []

        for E in total_energy:
            if E >= self.critical_thresholds['phase_4']:
                phases.append(4)
            elif E >= self.critical_thresholds['phase_3']:
                phases.append(3)
            elif E >= self.critical_thresholds['phase_2']:
                phases.append(2)
            else:
                phases.append(1)

        return total_energy, phases
