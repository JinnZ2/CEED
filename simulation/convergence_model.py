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
