# tests/test_convergence_model.py

import pytest
import numpy as np
from simulation.convergence_model import ConvergencePredictor

def test_model_initialization():
    predictor = ConvergencePredictor()
    assert isinstance(predictor.E_current, dict)
    assert all(k in predictor.E_current for k in ['solar', 'magnetic', 'atmospheric', 'oceanic'])

def test_energy_prediction_shape():
    predictor = ConvergencePredictor()
    t, solution = predictor.predict_convergence(years=1)
    assert len(t) == solution.shape[0]
    assert solution.shape[1] == 4  # 4 systems tracked

def test_phase_classification():
    predictor = ConvergencePredictor()
    t, solution = predictor.predict_convergence(years=1)
    total, phases = predictor.classify_phases(solution)
    assert len(phases) == len(total)
    assert all(1 <= p <= 4 for p in phases)

def test_total_energy_increase():
    predictor = ConvergencePredictor()
    t, solution = predictor.predict_convergence(years=1)
    total_energy = np.sum(solution, axis=1)
    assert total_energy[-1] > total_energy[0]  # If the world is heating up, this better be true
