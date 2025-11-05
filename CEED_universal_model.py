“””
Universal CEED Framework
Abstract feedback architecture applicable across domains

Domains: Climate, Finance, Ecosystems, Societies, AI Systems
“””

from dataclasses import dataclass
from typing import Callable, List, Tuple
import numpy as np

@dataclass
class FeedbackLoop:
“”“Represents a single feedback mechanism”””
name: str
polarity: str  # ‘positive’ or ‘negative’
strength: float
saturation_threshold: float = None

```
def apply(self, state: float, external_forcing: float = 0.0) -> float:
    """Apply feedback to current state"""
    if self.saturation_threshold and abs(state) > self.saturation_threshold:
        # Feedback saturates
        effective_strength = self.strength * (1 - (abs(state) / (self.saturation_threshold * 2)))
        effective_strength = max(0.1 * self.strength, effective_strength)
    else:
        effective_strength = self.strength
    
    if self.polarity == 'positive':
        return state * (1 + effective_strength) + external_forcing
    else:  # negative feedback
        return state * (1 - effective_strength) + external_forcing
```

@dataclass
class SystemState:
“”“Current state of a CEED system”””
energy: float  # Total system energy/stress level
retention: float  # How well system holds energy
dissipation: float  # How fast energy leaves
buffer_capacity: float  # Shock absorption remaining

```
def stability_metric(self) -> float:
    """
    Returns stability score
    > 1.0: System accumulating energy (unstable)
    ~ 1.0: Balanced
    < 1.0: System dissipating (stable)
    """
    return self.retention / self.dissipation if self.dissipation > 0 else float('inf')

def distance_to_tipping_point(self, threshold: float) -> float:
    """How close to critical threshold"""
    return (threshold - self.energy) / threshold
```

class CEEDSystem:
“””
Universal framework for modeling feedback-driven systems

```
Core Principle: Runaway occurs when:
1. Retention > Dissipation
2. Negative feedbacks saturate
3. Buffer capacity depletes
"""

def __init__(self, name: str):
    self.name = name
    self.feedbacks: List[FeedbackLoop] = []
    self.state = SystemState(
        energy=0.0,
        retention=1.0,
        dissipation=1.0,
        buffer_capacity=1.0
    )
    
    # Thresholds
    self.warning_threshold = 100.0
    self.critical_threshold = 200.0
    self.tipping_point = 300.0
    
def add_feedback(self, feedback: FeedbackLoop):
    """Add a feedback loop to the system"""
    self.feedbacks.append(feedback)

def compute_retention(self, energy: float) -> float:
    """
    Compute retention factor based on current energy
    
    Key insight: Retention typically drops at high energy
    (High-energy states are harder to maintain)
    """
    base_retention = 1.05  # Slight accumulation at low energy
    collapse_rate = 0.001  # How fast retention drops
    
    return base_retention * np.exp(-collapse_rate * energy**2)

def compute_dissipation(self, energy: float) -> float:
    """
    Compute dissipation based on current energy
    
    Typically increases nonlinearly with energy
    """
    linear_term = 0.05 * energy
    nonlinear_term = 0.001 * energy**1.5
    
    return linear_term + nonlinear_term

def apply_feedbacks(self, external_forcing: float = 0.0) -> float:
    """Apply all feedback loops to current state"""
    net_effect = self.state.energy
    
    for feedback in self.feedbacks:
        net_effect = feedback.apply(net_effect, external_forcing)
    
    return net_effect

def update(self, external_forcing: float = 0.0, dt: float = 0.1):
    """
    Update system state by one timestep
    
    Args:
        external_forcing: External shock/input
        dt: Timestep size
    """
    # Update retention and dissipation based on current energy
    self.state.retention = self.compute_retention(self.state.energy)
    self.state.dissipation = self.compute_dissipation(self.state.energy)
    
    # Apply feedbacks
    feedback_effect = self.apply_feedbacks(external_forcing)
    
    # Update energy
    energy_in = feedback_effect * self.state.retention
    energy_out = self.state.dissipation
    
    dE = (energy_in - energy_out) * dt
    self.state.energy += dE
    
    # Buffer capacity depletes with high energy states
    if self.state.energy > self.warning_threshold:
        depletion_rate = 0.01 * (self.state.energy / self.warning_threshold)
        self.state.buffer_capacity *= (1 - depletion_rate * dt)
        self.state.buffer_capacity = max(0.0, self.state.buffer_capacity)

def classify_state(self) -> str:
    """
    Classify current system state
    
    Returns:
        'stable', 'stressed', 'critical', or 'tipping'
    """
    if self.state.energy < self.warning_threshold:
        return 'stable'
    elif self.state.energy < self.critical_threshold:
        return 'stressed'
    elif self.state.energy < self.tipping_point:
        return 'critical'
    else:
        return 'tipping'

def diagnose(self) -> dict:
    """
    Generate diagnostic report
    
    Returns:
        Dictionary with system health metrics
    """
    return {
        'name': self.name,
        'state': self.classify_state(),
        'energy': self.state.energy,
        'stability': self.state.stability_metric(),
        'buffer_remaining': self.state.buffer_capacity,
        'distance_to_tipping': self.state.distance_to_tipping_point(self.tipping_point),
        'retention': self.state.retention,
        'dissipation': self.state.dissipation,
        'num_feedbacks': len(self.feedbacks)
    }

def simulate(self, steps: int, external_forcing_fn: Callable = None) -> List[dict]:
    """
    Run simulation for N steps
    
    Args:
        steps: Number of timesteps
        external_forcing_fn: Optional function(t) -> forcing
        
    Returns:
        List of diagnostic snapshots over time
    """
    history = []
    
    for t in range(steps):
        # Get external forcing if provided
        forcing = external_forcing_fn(t) if external_forcing_fn else 0.0
        
        # Update system
        self.update(external_forcing=forcing)
        
        # Record state
        snapshot = self.diagnose()
        snapshot['timestep'] = t
        history.append(snapshot)
    
    return history
```

# Example: Climate System

def create_climate_system() -> CEEDSystem:
“”“Create a CEED model of climate system”””
climate = CEEDSystem(“Climate”)

```
# Positive feedbacks
climate.add_feedback(FeedbackLoop(
    name="Water vapor",
    polarity="positive",
    strength=0.4,
    saturation_threshold=5.0  # Saturates at high temps
))

climate.add_feedback(FeedbackLoop(
    name="Ice-albedo",
    polarity="positive",
    strength=0.3,
    saturation_threshold=3.0  # Limited ice left to melt
))

# Negative feedbacks
climate.add_feedback(FeedbackLoop(
    name="Radiative cooling",
    polarity="negative",
    strength=0.5
))

climate.add_feedback(FeedbackLoop(
    name="Carbon sinks",
    polarity="negative",
    strength=0.3,
    saturation_threshold=4.0  # Sinks weaken with warming
))

return climate
```

# Example: Financial System

def create_financial_system() -> CEEDSystem:
“”“Create a CEED model of financial system”””
finance = CEEDSystem(“Finance”)

```
# Positive feedbacks
finance.add_feedback(FeedbackLoop(
    name="Leverage spiral",
    polarity="positive",
    strength=0.5,
    saturation_threshold=10.0
))

finance.add_feedback(FeedbackLoop(
    name="Panic selling",
    polarity="positive",
    strength=0.6,
    saturation_threshold=8.0
))

# Negative feedbacks
finance.add_feedback(FeedbackLoop(
    name="Central bank intervention",
    polarity="negative",
    strength=0.4,
    saturation_threshold=15.0  # Limited intervention capacity
))

finance.add_feedback(FeedbackLoop(
    name="Market liquidity",
    polarity="negative",
    strength=0.3,
    saturation_threshold=12.0  # Dries up in crisis
))

return finance
```

if **name** == “**main**”:
print(“CEED Universal Framework”)
print(”=” * 60)

```
# Create and test climate system
climate = create_climate_system()

print(f"\nInitial state: {climate.classify_state()}")
print(f"Stability metric: {climate.state.stability_metric():.2f}")

# Simulate with constant forcing
def forcing(t):
    return 5.0  # Constant external input

history = climate.simulate(steps=100, external_forcing_fn=forcing)

# Print final state
final = history[-1]
print(f"\nFinal state after 100 steps:")
print(f"  Classification: {final['state']}")
print(f"  Energy: {final['energy']:.1f}")
print(f"  Stability: {final['stability']:.2f}")
print(f"  Buffer remaining: {final['buffer_remaining']:.2%}")
print(f"  Distance to tipping point: {final['distance_to_tipping']:.2%}")
```
