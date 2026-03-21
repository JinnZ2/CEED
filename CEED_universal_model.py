"""
Universal CEED Framework
Domain-agnostic feedback architecture for modelling energy accumulation
in coupled systems (climate, finance, ecosystems, etc.).

Core ODE for total system energy E:

    dE/dt = F_ext(t) + sum_k f_k(E) - D(E)

where:
    F_ext(t) : external forcing [energy/time]
    f_k(E)   : k-th feedback contribution [energy/time]
             positive feedback: f_k = +s_k * E * sigma(E, E_sat_k)
             negative feedback: f_k = -s_k * E * sigma(E, E_sat_k)
    D(E)     : dissipation [energy/time]
    sigma(E, E_sat) = 1 / (1 + (E / E_sat)^2)  [dimensionless saturation]

Buffer capacity B depletes when E > E_warn:

    dB/dt = -b * (E / E_warn) * B   when E > E_warn, else 0
"""

from dataclasses import dataclass
from typing import Callable, List, Tuple, Optional
import numpy as np


@dataclass
class FeedbackLoop:
    """A single feedback mechanism.

    Attributes:
        name: Human-readable label.
        polarity: 'positive' (amplifying) or 'negative' (damping).
        strength: Feedback gain coefficient s_k [1/time].
        saturation_energy: Energy scale where feedback saturates [energy].
            None means no saturation (linear feedback).
    """
    name: str
    polarity: str
    strength: float
    saturation_energy: Optional[float] = None

    def rate(self, E: float) -> float:
        """Compute feedback contribution to dE/dt [energy/time].

        f_k(E) = (+/-) s_k * E * sigma(E)

        where sigma = 1/(1 + (E/E_sat)^2) if saturating, else 1.
        """
        sigma = 1.0
        if self.saturation_energy is not None and self.saturation_energy > 0:
            sigma = 1.0 / (1.0 + (E / self.saturation_energy) ** 2)

        sign = 1.0 if self.polarity == 'positive' else -1.0
        return sign * self.strength * E * sigma


@dataclass
class SystemState:
    """Instantaneous state of a CEED system.

    Attributes:
        energy: Total system energy [energy units].
        buffer_capacity: Remaining shock absorption capacity [0, 1].
    """
    energy: float
    buffer_capacity: float = 1.0


class CEEDSystem:
    """Universal framework for feedback-driven energy accumulation.

    The system evolves via:
        dE/dt = F_ext(t) + sum_k f_k(E) - D(E)

    Runaway occurs when positive feedbacks exceed dissipation:
        sum(positive f_k) > D(E) + sum(|negative f_k|)
    """

    def __init__(self, name: str):
        self.name = name
        self.feedbacks: List[FeedbackLoop] = []
        self.state = SystemState(energy=0.0, buffer_capacity=1.0)

        # Phase thresholds [energy units]
        self.warning_threshold = 100.0
        self.critical_threshold = 200.0
        self.tipping_point = 300.0

        # Dissipation parameters
        self.dissipation_linear = 0.05    # [1/time]
        self.dissipation_nonlinear = 0.001  # [1/(energy*time)]
        self.dissipation_exponent = 1.5   # nonlinear scaling power

        # Buffer depletion rate
        self.buffer_depletion_rate = 0.01  # [1/time]

    def add_feedback(self, feedback: FeedbackLoop):
        self.feedbacks.append(feedback)

    def dissipation_rate(self, E: float) -> float:
        """Total dissipation D(E) [energy/time].

        D(E) = alpha * E + beta * E^p

        where alpha is linear cooling, beta * E^p captures nonlinear
        radiative losses (motivated by Stefan-Boltzmann T^4 scaling).
        """
        linear = self.dissipation_linear * E
        nonlinear = self.dissipation_nonlinear * abs(E) ** self.dissipation_exponent
        return linear + nonlinear

    def total_feedback_rate(self, E: float) -> float:
        """Sum of all feedback contributions [energy/time]."""
        return sum(fb.rate(E) for fb in self.feedbacks)

    def dE_dt(self, E: float, external_forcing: float = 0.0) -> float:
        """Energy rate of change [energy/time].

        dE/dt = F_ext + sum_k f_k(E) - D(E)
        """
        return external_forcing + self.total_feedback_rate(E) - self.dissipation_rate(E)

    def update(self, external_forcing: float = 0.0, dt: float = 0.1):
        """Advance system by one timestep using forward Euler."""
        E = self.state.energy
        self.state.energy += self.dE_dt(E, external_forcing) * dt

        # Buffer depletion above warning threshold
        if self.state.energy > self.warning_threshold:
            rate = self.buffer_depletion_rate * (self.state.energy / self.warning_threshold)
            self.state.buffer_capacity *= (1.0 - rate * dt)
            self.state.buffer_capacity = max(0.0, self.state.buffer_capacity)

    def stability_metric(self) -> float:
        """Ratio of net amplification to dissipation.

        > 1.0: energy accumulating (unstable)
        ~ 1.0: balanced
        < 1.0: energy dissipating (stable)

        Returns inf if dissipation is zero.
        """
        E = self.state.energy
        D = self.dissipation_rate(E)
        if D == 0:
            return float('inf')
        net_feedback = self.total_feedback_rate(E)
        return (net_feedback + D) / D  # = 1 + feedback/dissipation

    def classify_state(self) -> str:
        """Classify current system state by energy level."""
        E = self.state.energy
        if E < self.warning_threshold:
            return 'stable'
        elif E < self.critical_threshold:
            return 'stressed'
        elif E < self.tipping_point:
            return 'critical'
        else:
            return 'tipping'

    def diagnose(self) -> dict:
        """Diagnostic snapshot of current system health."""
        E = self.state.energy
        return {
            'name': self.name,
            'state': self.classify_state(),
            'energy': E,
            'stability': self.stability_metric(),
            'buffer_remaining': self.state.buffer_capacity,
            'distance_to_tipping': (self.tipping_point - E) / self.tipping_point,
            'feedback_rate': self.total_feedback_rate(E),
            'dissipation_rate': self.dissipation_rate(E),
            'num_feedbacks': len(self.feedbacks),
        }

    def simulate(self, steps: int, dt: float = 0.1,
                 external_forcing_fn: Callable = None) -> List[dict]:
        """Run simulation for N steps.

        Args:
            steps: Number of timesteps.
            dt: Timestep size [time units].
            external_forcing_fn: Optional function(t) -> forcing [energy/time].

        Returns:
            List of diagnostic snapshots, one per timestep.
        """
        history = []
        for step in range(steps):
            t = step * dt
            forcing = external_forcing_fn(t) if external_forcing_fn else 0.0
            self.update(external_forcing=forcing, dt=dt)
            snapshot = self.diagnose()
            snapshot['timestep'] = step
            snapshot['time'] = t
            history.append(snapshot)
        return history


# ── Example systems ──────────────────────────────────────────────────

def create_climate_system() -> CEEDSystem:
    """Climate system with IPCC-motivated feedback strengths."""
    climate = CEEDSystem("Climate")

    # Positive feedbacks (amplify warming)
    climate.add_feedback(FeedbackLoop(
        name="Water vapor",
        polarity="positive",
        strength=0.04,          # [1/time]
        saturation_energy=500.0,
    ))
    climate.add_feedback(FeedbackLoop(
        name="Ice-albedo",
        polarity="positive",
        strength=0.03,
        saturation_energy=300.0,  # saturates as ice depletes
    ))

    # Negative feedbacks (stabilise)
    climate.add_feedback(FeedbackLoop(
        name="Radiative cooling",
        polarity="negative",
        strength=0.05,
    ))
    climate.add_feedback(FeedbackLoop(
        name="Carbon sinks",
        polarity="negative",
        strength=0.03,
        saturation_energy=400.0,  # sinks weaken with warming
    ))

    return climate


def create_financial_system() -> CEEDSystem:
    """Financial stress system with leverage and liquidity feedbacks."""
    finance = CEEDSystem("Finance")

    finance.add_feedback(FeedbackLoop(
        name="Leverage spiral",
        polarity="positive",
        strength=0.05,
        saturation_energy=1000.0,
    ))
    finance.add_feedback(FeedbackLoop(
        name="Panic selling",
        polarity="positive",
        strength=0.06,
        saturation_energy=800.0,
    ))
    finance.add_feedback(FeedbackLoop(
        name="Central bank intervention",
        polarity="negative",
        strength=0.04,
        saturation_energy=1500.0,
    ))
    finance.add_feedback(FeedbackLoop(
        name="Market liquidity",
        polarity="negative",
        strength=0.03,
        saturation_energy=1200.0,
    ))

    return finance


if __name__ == "__main__":
    print("CEED Universal Framework")
    print("=" * 60)

    climate = create_climate_system()
    climate.state.energy = 50.0  # start with moderate energy

    print(f"\nInitial state: {climate.classify_state()}")
    print(f"Stability metric: {climate.stability_metric():.3f}")

    def forcing(t):
        return 5.0  # constant external input [energy/time]

    history = climate.simulate(steps=100, dt=0.1, external_forcing_fn=forcing)

    final = history[-1]
    print(f"\nAfter 100 steps (t={final['time']:.1f}):")
    print(f"  Classification: {final['state']}")
    print(f"  Energy: {final['energy']:.1f}")
    print(f"  Stability metric: {final['stability']:.3f}")
    print(f"  Buffer remaining: {final['buffer_remaining']:.2%}")
    print(f"  Distance to tipping: {final['distance_to_tipping']:.2%}")
