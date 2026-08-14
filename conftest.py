"""
Pytest configuration.

Its presence at the repository root puts the root on sys.path during test
collection, so tests can do `from simulation.convergence_model import ...`
regardless of the directory pytest is invoked from.
"""
