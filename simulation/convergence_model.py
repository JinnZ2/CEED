"""
CEED Convergence Model
Multi-system energy balance with cross-domain coupling and anthropogenic forcing.

Each subsystem tracks a stored energy index E_i(t) governed by:

    dE_i/dt = S_i(t) + A_i(t)
              + (alpha_i(E,t) - lambda_i) * E_i
              - gamma_i * E_i^2
              + sum_j [eta_ij(E) * c_ij(E) * E_j]   (received, converted)
              - sum_j [c_ji(E) * E_i]                (sent, full amount)

where all parameters can be state-dependent:

    S_i(t)       : natural source/input rate [energy/yr]
    A_i(t)       : anthropogenic forcing — geological energy release [energy/yr]
    alpha_i(E,t) : effective retention rate [1/yr], increases with atmospheric
                   energy (greenhouse amplification) and anthropogenic load
    lambda_i     : radiative/dissipative loss rate [1/yr]
    gamma_i      : nonlinear dissipation [1/(energy*yr)] — 2nd law bound
    c_ij(E)      : transfer rate from j to i [1/yr], strengthens with
                   energy gradients (Clausius-Clapeyron, convective vigor)
    eta_ij(E)    : conversion efficiency [0-1], shifts with system state

Anthropogenic forcing:

    Fossil fuels represent ~300 Myr of stored photosynthetic energy released
    over ~200 yr — a 1.5-million-fold temporal compression.  This doesn't
    just add a source term; it modifies the system's own parameters:

    1. More CO2 → higher atmospheric retention (alpha_atm increases)
    2. Warmer atmosphere → stronger air-sea coupling (Clausius-Clapeyron)
    3. Warmer ocean → weaker carbon sinks → more CO2 → higher retention
    4. Rate of release is itself accelerating (~2%/yr historical growth)

    The anthropogenic load A(t) is split across atmospheric and oceanic
    subsystems (the two primary sinks for fossil energy).

Cross-system coupling is conservative (see previous docstring for details).
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
class AnthropogenicForcing:
    """Geological energy release from fossil fuel combustion.

    Models the temporal compression of ~300 Myr of stored solar energy
    into a ~200 yr release window, with accelerating extraction rate.

    The forcing modifies not just source terms but system parameters:
    retention rates, coupling strengths, and conversion efficiencies
    all shift as the anthropogenic energy load grows.

    Attributes:
        A_0: Current annual energy release rate [energy/yr].
            Calibrated to ~580 EJ/yr primary energy (2024), normalised
            to model units.
        growth_rate: Exponential growth rate of extraction [1/yr].
            Historical average ~2%/yr since 1950.
        peak_year: Year (from t=0) at which extraction peaks under
            resource constraints.  None = no peak (pure exponential).
        atm_fraction: Fraction of released energy deposited into
            atmospheric subsystem (rest goes to oceanic).
        alpha_sensitivity: How strongly atmospheric retention responds
            to anthropogenic load [dimensionless].
            delta_alpha = alpha_sensitivity * ln(1 + A(t)/A_0)
        coupling_sensitivity: How strongly coupling rates respond to
            energy gradients [dimensionless].
        efficiency_shift: How conversion efficiencies shift with
            total system energy [1/energy].
    """
    A_0: float = 8.0                # baseline release [energy/yr]
    growth_rate: float = 0.02       # 2%/yr exponential growth
    peak_year: Optional[float] = 50.0  # resource peak at t=50yr
    atm_fraction: float = 0.55      # 55% to atmosphere, 45% to ocean
    alpha_sensitivity: float = 0.02  # retention response to load
    coupling_sensitivity: float = 0.15  # coupling response to gradients
    efficiency_shift: float = 0.0003  # efficiency response to energy

    def release_rate(self, t: float) -> float:
        """Total anthropogenic energy release rate A(t) [energy/yr].

        Exponential growth up to peak_year, then logistic plateau:
            A(t) = A_0 * exp(r*t)                          if no peak
            A(t) = A_max / (1 + ((A_max/A_0) - 1)*exp(-r*t))  with peak

        The logistic form ensures smooth transition from exponential
        growth to resource-constrained plateau.
        """
        if self.peak_year is None:
            return self.A_0 * np.exp(self.growth_rate * t)

        # Logistic growth toward resource-constrained peak
        A_max = self.A_0 * np.exp(self.growth_rate * self.peak_year)
        ratio = A_max / self.A_0
        return A_max / (1.0 + (ratio - 1.0) * np.exp(-self.growth_rate * t))

    def cumulative_release(self, t: float) -> float:
        """Approximate cumulative energy released from t=0 to t [energy].

        For context: 300 Myr of storage released in ~200 yr means the
        cumulative release is a tiny fraction of geological storage but
        a massive perturbation to the current system.
        """
        if self.peak_year is None:
            return (self.A_0 / self.growth_rate) * (
                np.exp(self.growth_rate * t) - 1.0)

        # Approximate integral of logistic (exact for the form used)
        A_max = self.A_0 * np.exp(self.growth_rate * self.peak_year)
        ratio = A_max / self.A_0
        return (A_max / self.growth_rate) * np.log(
            (1.0 + (ratio - 1.0)) /
            (1.0 + (ratio - 1.0) * np.exp(-self.growth_rate * t))
        )


@dataclass
class SystemParameters:
    """Parameters for the coupled energy-balance ODE.

    All rates are per year.  Energy units are normalised (see baseline values).

    When anthropogenic forcing is enabled, retention rates, coupling
    strengths, and conversion efficiencies become state-dependent
    (functions of E and t), computed at each ODE evaluation.
    """

    # Initial energies [normalised energy units]
    E_initial: dict = field(default_factory=lambda: {
        'solar': 180.0,       # proxy: F10.7 index (sfu)
        'magnetic': 92.5,     # proxy: Kp-derived index
        'atmospheric': 118.0, # proxy: thermospheric density index
        'oceanic': 110.0,     # proxy: ocean heat content index
    })

    # Baseline retention rates [1/yr]
    # These are the pre-industrial / unperturbed values.
    # With anthropogenic forcing, effective alpha increases.
    alpha_retention: dict = field(default_factory=lambda: {
        'solar': 0.06,        # coronal magnetic confinement
        'magnetic': 0.025,    # ring current gradient drift
        'atmospheric': 0.09,  # greenhouse trapping (pre-industrial)
        'oceanic': 0.015,     # thermal inertia
    })

    # Linear dissipation rates [1/yr]
    lambda_dissipation: dict = field(default_factory=lambda: {
        'solar': 0.05,
        'magnetic': 0.02,
        'atmospheric': 0.08,
        'oceanic': 0.01,
    })

    # Nonlinear (quadratic) dissipation [1/(energy*yr)]
    gamma_nonlinear: dict = field(default_factory=lambda: {
        'solar': 0.0002,
        'magnetic': 0.0002,
        'atmospheric': 0.0002,
        'oceanic': 0.0002,
    })

    # Cross-system coupling: baseline values (transfer_rate, efficiency).
    # With anthropogenic forcing, both can be state-dependent.
    # Key (i, j) means "from j to i".
    coupling: dict = field(default_factory=lambda: {
        # Primary couplings
        ('magnetic', 'solar'):       (0.005, 0.15),
        ('atmospheric', 'solar'):    (0.003, 0.30),
        ('oceanic', 'solar'):        (0.001, 0.40),
        ('atmospheric', 'magnetic'): (0.002, 0.25),
        ('oceanic', 'atmospheric'):  (0.002, 0.35),
        # Secondary couplings (feedback)
        ('atmospheric', 'oceanic'):  (0.003, 0.30),
        ('magnetic', 'atmospheric'): (0.0005, 0.10),
        ('solar', 'magnetic'):       (0.0005, 0.10),
        ('magnetic', 'oceanic'):     (0.0001, 0.05),
    })

    # Phase classification thresholds on total energy [energy units]
    phase_thresholds: dict = field(default_factory=lambda: {
        'phase_1': 120,  # system stress
        'phase_2': 150,  # cross-system coupling onset
        'phase_3': 200,  # nonlinear amplification
        'phase_4': 300,  # cascade / collapse
    })

    # Anthropogenic forcing (None = disabled, static parameters)
    anthropogenic: Optional[AnthropogenicForcing] = None


class ConvergencePredictor:
    """Convergence model with optional anthropogenic forcing.

    When params.anthropogenic is None, all parameters are static
    (backward-compatible with previous model versions).

    When enabled, anthropogenic forcing:
    1. Adds a source term A_i(t) to atmospheric and oceanic subsystems
    2. Increases atmospheric retention (more GHGs → more trapping)
    3. Strengthens coupling between energetic subsystems
    4. Shifts conversion efficiencies with system state
    """

    def __init__(self, params: Optional[SystemParameters] = None):
        self.params = params or SystemParameters()
        self._sys_idx = {name: i for i, name in enumerate(SYSTEMS)}

    # ── Source terms ──────────────────────────────────────────────

    def source_rate(self, system: str, t: float) -> float:
        """Natural source/input rate S_i(t) [energy/yr].

        All sources are bounded to prevent unphysical growth.
        Secular trends use tanh saturation.
        """
        if system == 'solar':
            # 11-year sunspot cycle (periodic, naturally bounded)
            return 5.0 * (1.0 + 0.3 * np.cos(2 * np.pi * t / 11.0))
        elif system == 'magnetic':
            # Geomagnetic activity is driven by solar wind, not independent.
            # Base rate is small positive (inner core dynamo sustains field).
            # The solar coupling handles the solar-driven component.
            return 1.0 * (1.0 - 0.3 * np.tanh(0.01 * t))
        elif system == 'atmospheric':
            # Natural forcing trend, saturates
            return 3.0 * (1.0 + 0.5 * np.tanh(0.03 * t))
        else:  # oceanic
            # Ocean heat uptake trend, saturates
            return 1.0 * (1.0 + 0.3 * np.tanh(0.02 * t))

    def anthropogenic_source(self, system: str, t: float) -> float:
        """Anthropogenic energy injection A_i(t) [energy/yr].

        Geological stored energy is released primarily into the
        atmospheric and oceanic subsystems.  The atmospheric fraction
        represents direct radiative forcing from GHGs; the oceanic
        fraction represents ocean heat uptake.
        """
        anthro = self.params.anthropogenic
        if anthro is None:
            return 0.0

        A_total = anthro.release_rate(t)

        if system == 'atmospheric':
            return anthro.atm_fraction * A_total
        elif system == 'oceanic':
            return (1.0 - anthro.atm_fraction) * A_total
        else:
            return 0.0

    # ── State-dependent parameters ────────────────────────────────

    def effective_alpha(self, system: str, E: list, t: float) -> float:
        """Effective retention rate alpha_i(E, t) [1/yr].

        Anthropogenic forcing increases atmospheric retention:
            alpha_atm_eff = alpha_base * (1 + k * ln(1 + A(t)/A_0))

        Physical basis: more CO2 → more greenhouse trapping → more
        outgoing longwave radiation is absorbed and re-emitted downward.
        The logarithmic form matches the radiative forcing relationship.

        Oceanic retention also increases (thermal stratification from
        surface warming reduces vertical mixing, trapping heat):
            alpha_ocean_eff = alpha_base * (1 + k/2 * (E_ocean - E_ref)/E_ref)
        """
        p = self.params
        alpha_base = p.alpha_retention[system]
        anthro = p.anthropogenic

        if anthro is None:
            return alpha_base

        A_t = anthro.release_rate(t)

        if system == 'atmospheric':
            # Logarithmic response (mirrors CO2 radiative forcing)
            k = anthro.alpha_sensitivity
            return alpha_base * (1.0 + k * np.log(1.0 + A_t / anthro.A_0))

        elif system == 'oceanic':
            # Stratification feedback: warmer surface → less mixing → more retention
            E_ocean = E[self._sys_idx['oceanic']]
            E_ref = p.E_initial['oceanic']
            k = anthro.alpha_sensitivity * 0.5
            return alpha_base * (1.0 + k * max(0.0, E_ocean - E_ref) / E_ref)

        return alpha_base

    def effective_coupling(self, sys_i: str, sys_j: str,
                           E: list, t: float) -> Tuple[float, float]:
        """Effective coupling (rate, efficiency) for transfer j → i.

        Coupling strengthens with energy gradients between systems:
            c_eff = c_base * (1 + beta * |E_j - E_i| / E_ref)

        Physical basis: larger temperature/energy gradients drive
        stronger fluxes (Clausius-Clapeyron for moisture, Fourier's
        law for heat, Ohm's law for currents).

        Conversion efficiency shifts with total system energy:
            eta_eff = eta_base * (1 + delta * (E_total - E_ref) / E_ref)
            clamped to [0.01, 0.95] to respect thermodynamic limits.

        At higher total energy, conversion processes are more vigorous
        (higher T → faster reaction kinetics, more turbulent mixing)
        but capped below 1.0 (no perpetual motion).
        """
        p = self.params
        entry = p.coupling.get((sys_i, sys_j))
        if entry is None:
            return (0.0, 0.0)

        c_base, eta_base = entry
        anthro = p.anthropogenic

        if anthro is None:
            return (c_base, eta_base)

        i = self._sys_idx[sys_i]
        j = self._sys_idx[sys_j]
        E_total = sum(E)
        E_ref = sum(p.E_initial.values())

        # Coupling rate scales with gradient
        beta = anthro.coupling_sensitivity
        gradient = abs(E[j] - E[i])
        E_scale = max(p.E_initial[sys_j], 1.0)
        c_eff = c_base * (1.0 + beta * gradient / E_scale)

        # Efficiency shifts with total energy
        delta = anthro.efficiency_shift
        eta_eff = eta_base * (1.0 + delta * (E_total - E_ref))
        eta_eff = max(0.01, min(0.95, eta_eff))

        return (c_eff, eta_eff)

    # ── ODE ───────────────────────────────────────────────────────

    def energy_derivative(self, E: list, t: float) -> list:
        """RHS of the coupled ODE system.

        dE_i/dt = S_i(t) + A_i(t)
                  + (alpha_i(E,t) - lambda_i) * E_i
                  - gamma_i * E_i^2
                  + sum_j [eta_ij(E)*c_ij(E)*E_j]   (received)
                  - sum_j [c_ji(E)*E_i]               (sent)
        """
        p = self.params
        dE_dt = []

        for i, sys_i in enumerate(SYSTEMS):
            E_i = E[i]
            alpha = self.effective_alpha(sys_i, E, t)
            lam = p.lambda_dissipation[sys_i]
            gam = p.gamma_nonlinear[sys_i]

            source = self.source_rate(sys_i, t)
            anthro = self.anthropogenic_source(sys_i, t)

            # Retention and dissipation only act on positive energy.
            # At E=0 the system is empty — nothing to retain or dissipate.
            # Soft floor via max(E_i, 0) prevents negative-energy runaway
            # where gamma*E^2 would push further negative.
            E_pos = max(E_i, 0.0)
            net_retention = (alpha - lam) * E_pos
            nonlinear_loss = gam * E_pos ** 2

            # If E_i went negative (unphysical), add restoring force
            # to push back toward zero.  This is a soft floor.
            restoring = 0.0
            if E_i < 0:
                restoring = -E_i  # linear spring toward zero

            # Energy received from other systems (only from positive E)
            coupling_in = 0.0
            for j, sys_j in enumerate(SYSTEMS):
                if j != i:
                    c_ij, eta_ij = self.effective_coupling(sys_i, sys_j, E, t)
                    coupling_in += eta_ij * c_ij * max(E[j], 0.0)

            # Energy sent to other systems (only if we have energy)
            coupling_out = 0.0
            for j, sys_j in enumerate(SYSTEMS):
                if j != i:
                    c_ji, _ = self.effective_coupling(sys_j, sys_i, E, t)
                    coupling_out += c_ji * E_pos

            dE_dt.append(source + anthro + net_retention - nonlinear_loss
                         + coupling_in - coupling_out + restoring)

        return dE_dt

    # ── Integration ───────────────────────────────────────────────

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
        """Classify each timestep into a phase based on total energy."""
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
    remains deterministic.
    """

    def __init__(self, params: Optional[SystemParameters] = None):
        super().__init__(params)
        self.meteor_rate = 0.5
        self.meteor_energy_range = (1.0, 15.0)
        self.launches_per_year = 150
        self.launch_energy = 2.0
        self.sink_baseline = 0.15
        self.sink_saturation_energy = 800.0
        self.events: List[ExternalEvent] = []

    def _generate_events(self, years: float) -> List[ExternalEvent]:
        """Pre-generate all stochastic external events."""
        events = []

        n_meteors = np.random.poisson(self.meteor_rate * years)
        for _ in range(n_meteors):
            t = np.random.uniform(0, years)
            energy = np.random.uniform(*self.meteor_energy_range)
            events.append(ExternalEvent(t, energy, 'meteor'))

        n_launches = int(self.launches_per_year * years)
        for i in range(n_launches):
            t = (i / self.launches_per_year) + np.random.uniform(-0.01, 0.01)
            if 0 <= t < years:
                events.append(ExternalEvent(t, self.launch_energy, 'launch'))

        n_volcanic = np.random.poisson(0.6 * years)
        for _ in range(n_volcanic):
            t = np.random.uniform(0, years)
            energy = np.random.uniform(5, 25)
            events.append(ExternalEvent(t, -energy, 'volcanic'))

        return sorted(events, key=lambda e: e.time)

    def _event_forcing(self, t: float) -> float:
        """Sum external event contributions at time t (Gaussian pulses)."""
        total = 0.0
        for event in self.events:
            dt = (t - event.time) / event.duration
            if abs(dt) < 5:
                total += (event.energy / (event.duration * np.sqrt(2 * np.pi))
                          * np.exp(-0.5 * dt ** 2))
        return total

    def _unknown_sink_rate(self, E_total: float) -> float:
        """Saturating dissipation from uncharacterised geophysical sinks."""
        sat = max(0.1, 1.0 - E_total / self.sink_saturation_energy)
        return self.sink_baseline * E_total * sat

    def energy_derivative(self, E: list, t: float) -> list:
        """RHS with external forcing and unknown sinks added."""
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
        """Integrate with optional external events."""
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
                continue
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
    print("CEED Convergence Model")
    print("=" * 60)

    # ── Baseline (no anthropogenic forcing) ──
    params_base = SystemParameters()
    base = ConvergencePredictor(params_base)
    t_b, sol_b = base.predict_convergence(years=10)
    E_b, ph_b = base.classify_phases(sol_b)

    # ── With anthropogenic forcing ──
    params_anthro = SystemParameters(
        anthropogenic=AnthropogenicForcing()
    )
    anthro = ConvergencePredictor(params_anthro)
    t_a, sol_a = anthro.predict_convergence(years=10)
    E_a, ph_a = anthro.classify_phases(sol_a)

    print(f"\nBASELINE (no anthropogenic):")
    print(f"  E_total: {E_b[0]:.1f} -> {E_b[-1]:.1f}")
    print(f"  Peak phase: {max(ph_b)}")

    print(f"\nWITH ANTHROPOGENIC FORCING:")
    print(f"  E_total: {E_a[0]:.1f} -> {E_a[-1]:.1f}")
    print(f"  Peak phase: {max(ph_a)}")
    af = params_anthro.anthropogenic
    print(f"  A(t=0):  {af.release_rate(0):.2f} energy/yr")
    print(f"  A(t=10): {af.release_rate(10):.2f} energy/yr")
    print(f"  Cumulative release: {af.cumulative_release(10):.1f} energy")

    # Show how parameters shifted
    E_final = list(sol_a[-1])
    alpha_base = params_anthro.alpha_retention['atmospheric']
    alpha_eff = anthro.effective_alpha('atmospheric', E_final, 10.0)
    print(f"\n  Atmospheric retention: {alpha_base:.4f} -> {alpha_eff:.4f} "
          f"(+{(alpha_eff/alpha_base - 1)*100:.1f}%)")

    c_base, eta_base = params_anthro.coupling[('oceanic', 'atmospheric')]
    c_eff, eta_eff = anthro.effective_coupling('oceanic', 'atmospheric',
                                                E_final, 10.0)
    print(f"  Atm->Ocean coupling:  {c_base:.4f} -> {c_eff:.4f} "
          f"(+{(c_eff/c_base - 1)*100:.1f}%)")
    print(f"  Atm->Ocean efficiency: {eta_base:.2f} -> {eta_eff:.2f}")

    # Show per-subsystem energy at end
    print(f"\n  Final subsystem energies:")
    for i, s in enumerate(SYSTEMS):
        print(f"    {s:>12s}: {sol_b[-1, i]:7.1f} (base) -> "
              f"{sol_a[-1, i]:7.1f} (anthro)  "
              f"delta={sol_a[-1, i] - sol_b[-1, i]:+.1f}")
