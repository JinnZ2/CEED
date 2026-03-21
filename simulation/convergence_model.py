"""
CEED Convergence Model
Multi-system energy balance with cross-domain coupling.

Each subsystem tracks a stored energy index E_i(t) governed by:

    dE_i/dt = S_i(t) + alpha_i * E_i - lambda_i * E_i - gamma_i * E_i^2
              + sum_j(c_ij * E_j) + F_ext_i(t)

which can be written:

    dE_i/dt = S_i(t) + (alpha_i - lambda_i) * E_i - gamma_i * E_i^2
              + sum_j(c_ij * E_j) + F_ext_i(t)

where:
    S_i(t)       : external source/input rate [energy/yr]
    alpha_i      : retention rate — rate at which the system reinforces
                   its own energy state [1/yr].  Physical mechanisms:
                   greenhouse trapping (atm), magnetic confinement (mag),
                   ocean thermal inertia (oceanic), coronal storage (solar)
    lambda_i     : radiative/dissipative loss rate [1/yr]
    gamma_i      : nonlinear dissipation [1/(energy*yr)] — dominates at
                   high E, ensuring bounded solutions (2nd law: entropy
                   production increases with energy, enhanced radiation)
    c_ij         : cross-system coupling [1/yr]
    F_ext_i(t)   : external forcing (meteors, launches, volcanic) [energy/yr]

Energy conservation (1st law): energy entering the system either stays
(retention) or leaves (dissipation).  The net rate (alpha - lambda) controls
whether the system accumulates or loses energy.  CEED's thesis: when
alpha_i >= lambda_i across multiple coupled subsystems, energy accumulates
until nonlinear dissipation (gamma * E^2) restores balance — or the system
transitions to a new state.

2nd law constraint: gamma_i > 0 guarantees that dissipation always wins
at sufficiently high E (entropy production scales superlinearly with
energy flux), preventing unphysical unbounded growth.
"""

import numpy as np
from scipy.integrate import odeint
from dataclasses import dataclass, field
from typing import List, Tuple, Optional


SYSTEMS = ['solar', 'magnetic', 'atmospheric', 'oceanic']
N_SYSTEMS = len(SYSTEMS)


@dataclass
class ExternalEvent:
    """An impulsive external energy perturbation.

    Attributes:
        time: Event time [years].
        energy: Energy deposited (positive) or released (negative) [energy units].
        event_type: Category label for the event.
        duration: Characteristic duration [years]. Default 1/365 (~1 day).
    """
    time: float
    energy: float
    event_type: str
    duration: float = 1.0 / 365.0


@dataclass
class SystemParameters:
    """Parameters for the coupled energy-balance ODE.

    All rates are per year.  Energy units are normalised (see baseline values).
    """

    # Initial energies [normalised energy units]
    E_initial: dict = field(default_factory=lambda: {
        'solar': 180.0,       # proxy: F10.7 index (sfu)
        'magnetic': 92.5,     # proxy: Kp-derived index
        'atmospheric': 118.0, # proxy: thermospheric density index
        'oceanic': 110.0,     # proxy: ocean heat content index
    })

    # Retention rates [1/yr]
    # Rate at which each subsystem reinforces its own energy state.
    # Physical basis:
    #   solar       — coronal magnetic confinement of plasma
    #   magnetic    — ring current self-sustaining via gradient drift
    #   atmospheric — greenhouse trapping of outgoing longwave radiation
    #   oceanic     — thermal inertia of deep ocean heat storage
    alpha_retention: dict = field(default_factory=lambda: {
        'solar': 0.06,
        'magnetic': 0.025,
        'atmospheric': 0.09,
        'oceanic': 0.015,
    })

    # Linear dissipation rates [1/yr]
    # Radiative losses, particle precipitation, heat flux to space.
    # Net growth when alpha > lambda; net decay when alpha < lambda.
    lambda_dissipation: dict = field(default_factory=lambda: {
        'solar': 0.05,
        'magnetic': 0.02,
        'atmospheric': 0.08,
        'oceanic': 0.01,
    })

    # Nonlinear (quadratic) dissipation [1/(energy*yr)]
    # Enforces 2nd-law bound: at high E, entropy production and
    # radiative losses scale superlinearly, preventing runaway.
    gamma_nonlinear: dict = field(default_factory=lambda: {
        'solar': 0.0002,
        'magnetic': 0.0002,
        'atmospheric': 0.0002,
        'oceanic': 0.0002,
    })

    # Cross-system coupling matrix c_ij [1/yr]
    # c_ij > 0 means energy in system j drives growth in system i.
    # Only physically motivated couplings are nonzero.
    coupling: dict = field(default_factory=lambda: {
        # solar -> magnetic: geomagnetic storms driven by solar wind
        ('magnetic', 'solar'): 0.005,
        # solar -> atmospheric: EUV heating of thermosphere
        ('atmospheric', 'solar'): 0.003,
        # magnetic -> atmospheric: Joule heating from auroral currents
        ('atmospheric', 'magnetic'): 0.002,
        # atmospheric -> oceanic: air-sea heat flux
        ('oceanic', 'atmospheric'): 0.001,
    })

    # Phase classification thresholds on total energy [energy units]
    phase_thresholds: dict = field(default_factory=lambda: {
        'phase_1': 120,  # system stress
        'phase_2': 150,  # cross-system coupling onset
        'phase_3': 200,  # nonlinear amplification
        'phase_4': 300,  # cascade / collapse
    })


class ConvergencePredictor:
    """Baseline convergence model (no external events)."""

    def __init__(self, params: Optional[SystemParameters] = None):
        self.params = params or SystemParameters()
        self._sys_idx = {name: i for i, name in enumerate(SYSTEMS)}

    def source_rate(self, system: str, t: float) -> float:
        """External source/input rate S_i(t) [energy/yr].

        Solar: 11-year sunspot cycle modulation (F10.7 proxy).
        Magnetic: secular geomagnetic field weakening trend.
        Atmospheric: anthropogenic greenhouse gas trend.
        Oceanic: anthropogenic ocean heat uptake trend.
        """
        if system == 'solar':
            return 5.0 * (1.0 + 0.3 * np.cos(2 * np.pi * t / 11.0))
        elif system == 'magnetic':
            return -2.0 * (1.0 + 0.1 * t)
        elif system == 'atmospheric':
            return 3.0 * (1.0 + 0.05 * t)
        else:  # oceanic
            return 1.0 * (1.0 + 0.02 * t)

    def energy_derivative(self, E: list, t: float) -> list:
        """RHS of the coupled ODE system.

        dE_i/dt = S_i(t) + (alpha_i - lambda_i)*E_i - gamma_i*E_i^2
                  + sum_j c_ij*E_j

        The (alpha - lambda) term is the net retention: positive means the
        system accumulates energy faster than it radiates.  The gamma*E^2
        term ensures dissipation always wins at high E (2nd law).
        """
        p = self.params
        dE_dt = []

        for i, sys_i in enumerate(SYSTEMS):
            E_i = E[i]
            alpha = p.alpha_retention[sys_i]
            lam = p.lambda_dissipation[sys_i]
            gam = p.gamma_nonlinear[sys_i]

            source = self.source_rate(sys_i, t)
            net_retention = (alpha - lam) * E_i   # >0 when accumulating
            nonlinear_loss = gam * E_i ** 2        # always positive, grows with E

            # Cross-system coupling: sum over j != i
            coupling_gain = 0.0
            for j, sys_j in enumerate(SYSTEMS):
                if j != i:
                    c_ij = p.coupling.get((sys_i, sys_j), 0.0)
                    coupling_gain += c_ij * E[j]

            dE_dt.append(source + net_retention - nonlinear_loss + coupling_gain)

        return dE_dt

    def predict_convergence(self, years: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
        """Integrate the ODE forward in time.

        Returns:
            t: time array [years], shape (N,)
            solution: energy array, shape (N, 4)
        """
        n_steps = max(int(years * 12), 2)
        t = np.linspace(0, years, n_steps)
        E0 = [self.params.E_initial[s] for s in SYSTEMS]
        solution = odeint(self.energy_derivative, E0, t)
        return t, solution

    def classify_phases(self, solution: np.ndarray) -> Tuple[np.ndarray, list]:
        """Classify each timestep into a phase based on total energy.

        Returns:
            total_energy: shape (N,)
            phases: list of int (1-4)
        """
        thresholds = self.params.phase_thresholds
        total_energy = np.sum(solution, axis=1)
        phases = []
        for E in total_energy:
            if E >= thresholds['phase_4']:
                phases.append(4)
            elif E >= thresholds['phase_3']:
                phases.append(3)
            elif E >= thresholds['phase_2']:
                phases.append(2)
            else:
                phases.append(1)
        return total_energy, phases


class ExtendedConvergencePredictor(ConvergencePredictor):
    """Convergence model with external stochastic forcing.

    External events are pre-generated before integration so the ODE RHS
    remains deterministic (required for reliable ODE solver behaviour).

    Additional dissipation from uncharacterised sinks (volcanic, seismic)
    is modelled as a saturating linear term:

        D_unk(E_tot) = mu * E_tot * max(0.1, 1 - E_tot / E_sat)

    where mu is the baseline sink fraction and E_sat is the saturation energy.
    """

    def __init__(self, params: Optional[SystemParameters] = None):
        super().__init__(params)

        # Meteor impact parameters
        self.meteor_rate = 0.5            # events/yr (Poisson rate)
        self.meteor_energy_range = (1.0, 15.0)  # [energy units]

        # Satellite launch parameters
        self.launches_per_year = 150      # global launch cadence
        self.launch_energy = 2.0          # ionospheric perturbation per launch

        # Unknown sink parameters
        self.sink_baseline = 0.15         # fractional loss rate at low E [1/yr]
        self.sink_saturation_energy = 800.0  # energy where sinks saturate

        # Pre-generated events (populated by predict_convergence)
        self.events: List[ExternalEvent] = []

    def _generate_events(self, years: float) -> List[ExternalEvent]:
        """Pre-generate all stochastic external events."""
        events = []

        # Meteor / bolide impacts (Poisson process)
        n_meteors = np.random.poisson(self.meteor_rate * years)
        for _ in range(n_meteors):
            t = np.random.uniform(0, years)
            energy = np.random.uniform(*self.meteor_energy_range)
            events.append(ExternalEvent(t, energy, 'meteor'))

        # Satellite launches (quasi-periodic with jitter)
        n_launches = int(self.launches_per_year * years)
        for i in range(n_launches):
            t = (i / self.launches_per_year) + np.random.uniform(-0.01, 0.01)
            if 0 <= t < years:
                events.append(ExternalEvent(t, self.launch_energy, 'launch'))

        # Volcanic / seismic dissipation events (Poisson, ~5% chance per month)
        n_volcanic = np.random.poisson(0.6 * years)  # ~0.6 events/yr
        for _ in range(n_volcanic):
            t = np.random.uniform(0, years)
            energy = np.random.uniform(5, 25)
            events.append(ExternalEvent(t, -energy, 'volcanic'))

        return sorted(events, key=lambda e: e.time)

    def _event_forcing(self, t: float) -> float:
        """Sum external event contributions at time t.

        Each event is modelled as a Gaussian pulse with width = event.duration:
            F(t) = energy / (duration * sqrt(2*pi)) * exp(-0.5*((t - t0)/duration)^2)
        This avoids the original approach of point-matching with a tolerance,
        which missed events between ODE solver adaptive steps.
        """
        total = 0.0
        for event in self.events:
            dt = (t - event.time) / event.duration
            if abs(dt) < 5:  # truncate at 5 sigma for speed
                total += (event.energy / (event.duration * np.sqrt(2 * np.pi))
                          * np.exp(-0.5 * dt ** 2))
        return total

    def _unknown_sink_rate(self, E_total: float) -> float:
        """Saturating dissipation from uncharacterised geophysical sinks.

        Returns a dissipation rate [energy/yr] that is subtracted from dE/dt.
        """
        sat = max(0.1, 1.0 - E_total / self.sink_saturation_energy)
        return self.sink_baseline * E_total * sat

    def energy_derivative(self, E: list, t: float) -> list:
        """RHS with external forcing and unknown sinks added.

        When no events are loaded (include_external=False), this reduces
        exactly to the base ConvergencePredictor ODE.
        """
        dE_dt = super().energy_derivative(E, t)

        if not self.events:
            return dE_dt

        E_total = sum(E)
        ext = self._event_forcing(t)
        sink = self._unknown_sink_rate(E_total)

        per_system_ext = ext / N_SYSTEMS
        per_system_sink = sink / N_SYSTEMS

        for i in range(N_SYSTEMS):
            dE_dt[i] += per_system_ext - per_system_sink

        return dE_dt

    def predict_convergence(self, years: float = 3.0,
                            include_external: bool = True
                            ) -> Tuple[np.ndarray, np.ndarray, List[ExternalEvent]]:
        """Integrate with optional external events.

        Returns:
            t: time array [years]
            solution: energy array, shape (N, 4)
            events: list of ExternalEvent used in this run
        """
        if include_external:
            self.events = self._generate_events(years)
        else:
            self.events = []

        n_steps = max(int(years * 12), 2)
        t = np.linspace(0, years, n_steps)
        E0 = [self.params.E_initial[s] for s in SYSTEMS]
        solution = odeint(self.energy_derivative, E0, t)
        return t, solution, self.events

    def analyze_event_timing(self, t: np.ndarray, solution: np.ndarray,
                             events: List[ExternalEvent]) -> list:
        """Identify external events that coincide with high-energy states."""
        total_energy = np.sum(solution, axis=1)
        threshold = self.params.phase_thresholds['phase_2']

        critical = []
        for event in events:
            if event.energy <= 0:
                continue  # skip dissipation events
            idx = np.argmin(np.abs(t - event.time))
            E_at = total_energy[idx]
            if E_at > threshold:
                critical.append({
                    'event': event,
                    'system_energy': E_at,
                    'amplification_risk': E_at / threshold,
                })
        return critical


if __name__ == "__main__":
    predictor = ExtendedConvergencePredictor()

    print("CEED Convergence Model")
    print("=" * 50)

    # Baseline run (no external events)
    t_base, sol_base, _ = predictor.predict_convergence(years=3, include_external=False)
    E_base, phases_base = predictor.classify_phases(sol_base)

    # Run with external events
    t_ext, sol_ext, events = predictor.predict_convergence(years=3, include_external=True)
    E_ext, phases_ext = predictor.classify_phases(sol_ext)

    print(f"\nBASELINE (no external events):")
    print(f"  Start: {E_base[0]:.1f} energy units")
    print(f"  End:   {E_base[-1]:.1f} energy units")
    print(f"  Peak phase: {max(phases_base)}")

    print(f"\nWITH EXTERNAL EVENTS:")
    print(f"  Start: {E_ext[0]:.1f} energy units")
    print(f"  End:   {E_ext[-1]:.1f} energy units")
    print(f"  Peak phase: {max(phases_ext)}")
    print(f"  Total events: {len(events)}")

    critical = predictor.analyze_event_timing(t_ext, sol_ext, events)
    if critical:
        print(f"\nCRITICAL TIMING EVENTS: {len(critical)}")
        for c in critical[:5]:
            print(f"  {c['event'].event_type} at t={c['event'].time:.2f}yr: "
                  f"E_total={c['system_energy']:.1f} "
                  f"({c['amplification_risk']:.2f}x phase-2 threshold)")
