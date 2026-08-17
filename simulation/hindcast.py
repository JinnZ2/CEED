"""
Hindcast Validation Framework for the CEED Convergence Model.

Runs the model against historical observations and scores the results.
Any AI (or human) can use this to:

1. Validate model accuracy against real data
2. Tune parameters to improve fit
3. Compare model variants (baseline vs anthropogenic, different couplings)
4. Understand where the model succeeds and fails

Usage:
    python -m simulation.hindcast                    # run all hindcasts
    python -m simulation.hindcast --system solar     # single subsystem
    python -m simulation.hindcast --plot             # show plots

Scoring metrics:
    RMSE  : Root Mean Squared Error (lower is better)
    MAE   : Mean Absolute Error (lower is better)
    R²    : Coefficient of determination (1.0 = perfect, 0 = random)
    NRMSE : Normalised RMSE (fraction of observed range)
    Bias  : Mean signed error (positive = model runs hot)

A good hindcast should have:
    R² > 0.5   (explains more than half the variance)
    NRMSE < 0.3 (error < 30% of observed range)
    |Bias| small relative to RMSE
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import argparse

from simulation.convergence_model import (
    ConvergencePredictor,
    SystemParameters,
    AnthropogenicForcing,
    SYSTEMS,
)
from simulation.unit_bridge import (
    SCALES,
    ALL_OBSERVED,
    REFERENCE_YEAR,
    observations_to_energy,
    energy_to_observations,
)


@dataclass
class HindcastScore:
    """Validation metrics for a single subsystem hindcast.

    All metrics are computed in physical units (not energy indices)
    so they're interpretable.
    """
    system: str
    units: str
    rmse: float
    mae: float
    r_squared: float
    nrmse: float       # RMSE / range of observations
    bias: float         # mean(predicted - observed)
    n_points: int
    years: List[int]
    observed: np.ndarray
    predicted: np.ndarray

    def grade(self) -> str:
        """Letter grade based on R² and NRMSE."""
        if self.r_squared > 0.8 and self.nrmse < 0.15:
            return 'A'
        elif self.r_squared > 0.6 and self.nrmse < 0.25:
            return 'B'
        elif self.r_squared > 0.3 and self.nrmse < 0.4:
            return 'C'
        elif self.r_squared > 0.0:
            return 'D'
        else:
            return 'F'

    def summary(self) -> str:
        """One-line summary."""
        return (f"{self.system:<14s} R²={self.r_squared:+.3f}  "
                f"RMSE={self.rmse:.2f} {self.units:<10s} "
                f"NRMSE={self.nrmse:.2f}  Bias={self.bias:+.2f}  "
                f"Grade={self.grade()}")


def compute_scores(system: str, years: np.ndarray,
                   observed: np.ndarray,
                   predicted: np.ndarray) -> HindcastScore:
    """Compute validation metrics between observed and predicted.

    Args:
        system: Subsystem name.
        years: Calendar years.
        observed: Observed values in physical units.
        predicted: Predicted values in physical units.

    Returns:
        HindcastScore with all metrics.
    """
    units = SCALES[system].units
    residuals = predicted - observed
    n = len(observed)

    rmse = np.sqrt(np.mean(residuals ** 2))
    mae = np.mean(np.abs(residuals))
    bias = np.mean(residuals)

    obs_range = np.max(observed) - np.min(observed)
    nrmse = rmse / obs_range if obs_range > 0 else float('inf')

    obs_mean = np.mean(observed)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((observed - obs_mean) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return HindcastScore(
        system=system, units=units, rmse=rmse, mae=mae,
        r_squared=r_squared, nrmse=nrmse, bias=bias,
        n_points=n, years=list(years.astype(int)),
        observed=observed, predicted=predicted,
    )


def integrate_from(start_year: int,
                   years: float,
                   params: SystemParameters
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """Initialise from observations at start_year and integrate forward.

    Shared by the hindcast and the forecast so the two cannot drift apart:
    a forecast that initialises differently from the hindcast that validated
    it is not testing the same model.

    The model's internal clock has t=0 at REFERENCE_YEAR, so the time-varying
    terms are shifted by (start_year - REFERENCE_YEAR). `params` is mutated in
    place to carry the new initial conditions.

    Args:
        start_year: Calendar year to initialise from. Must have observations.
        years: How far forward to integrate.
        params: Model parameters. E_initial is overwritten.

    Returns:
        t (years since start_year), solution of shape (N, 4).
    """
    for sys_name in SYSTEMS:
        obs_data = ALL_OBSERVED[sys_name]
        if start_year in obs_data:
            params.E_initial[sys_name] = SCALES[sys_name].obs_to_energy(
                obs_data[start_year])

    predictor = ConvergencePredictor(params)
    offset = start_year - REFERENCE_YEAR

    original_source = predictor.source_rate
    original_anthro = predictor.anthropogenic_source
    original_alpha = predictor.effective_alpha

    predictor.source_rate = lambda system, t: original_source(system, t + offset)
    predictor.anthropogenic_source = (
        lambda system, t: original_anthro(system, t + offset))
    predictor.effective_alpha = (
        lambda system, E, t: original_alpha(system, E, t + offset))

    return predictor.predict_convergence(years=float(years))


def run_hindcast(params: Optional[SystemParameters] = None,
                 systems: Optional[List[str]] = None,
                 start_year: int = 2010,
                 end_year: int = 2024,
                 score_from: Optional[int] = None,
                 score_to: Optional[int] = None,
                 ) -> Dict[str, HindcastScore]:
    """Run hindcast validation for specified subsystems.

    The model is initialised at start_year using observed data for that year,
    then integrated forward.  Predictions are compared to observations at
    each subsequent year.

    score_from/score_to restrict SCORING to a sub-window without changing the
    integration, which is what makes a held-out test possible: the trajectory
    still starts at start_year and receives no observations after it, so
    scoring a later window is a genuine free-running forecast rather than a
    re-initialised one. See simulation/validation.py.

    Args:
        params: Model parameters (default: SystemParameters with anthropogenic).
        systems: List of subsystem names to validate. Default: all.
        start_year: First year of hindcast period (model initialised here).
        end_year: Last year integrated.
        score_from: First year to score. Default start_year + 1.
        score_to: Last year to score. Default end_year.

    Returns:
        Dict mapping system name to HindcastScore.
    """
    if systems is None:
        systems = list(SYSTEMS)

    if params is None:
        params = SystemParameters(anthropogenic=AnthropogenicForcing())

    t, solution = integrate_from(start_year, end_year - start_year, params)

    # Extract predictions at observation years.
    #
    # start_year is EXCLUDED from scoring: the model is initialised from the
    # observation at that year, so its residual there is identically zero by
    # construction. Including it donated one free perfect point per subsystem
    # and inflated R².
    scores = {}
    for sys_name in systems:
        obs_data = ALL_OBSERVED[sys_name]
        lo = start_year + 1 if score_from is None else max(score_from, start_year + 1)
        hi = end_year if score_to is None else min(score_to, end_year)
        obs_years = sorted([yr for yr in obs_data if lo <= yr <= hi])
        if len(obs_years) < 2:
            continue

        obs_vals = np.array([obs_data[yr] for yr in obs_years])
        sys_idx = SYSTEMS.index(sys_name)

        # Interpolate model solution at observation times
        t_obs = np.array([float(yr - start_year) for yr in obs_years])
        E_interp = np.interp(t_obs, t, solution[:, sys_idx])

        # Convert to physical units
        pred_obs = energy_to_observations(sys_name, E_interp)

        scores[sys_name] = compute_scores(
            sys_name,
            np.array(obs_years),
            obs_vals,
            pred_obs,
        )

    return scores


def run_comparison(start_year: int = 2010, end_year: int = 2024
                   ) -> Tuple[Dict[str, HindcastScore],
                              Dict[str, HindcastScore]]:
    """Run baseline vs anthropogenic hindcast comparison.

    Returns:
        (baseline_scores, anthropogenic_scores)
    """
    baseline_params = SystemParameters()
    anthro_params = SystemParameters(anthropogenic=AnthropogenicForcing())

    baseline_scores = run_hindcast(baseline_params,
                                   start_year=start_year,
                                   end_year=end_year)
    anthro_scores = run_hindcast(anthro_params,
                                 start_year=start_year,
                                 end_year=end_year)

    return baseline_scores, anthro_scores


def print_report(scores: Dict[str, HindcastScore],
                 title: str = "Hindcast Report"):
    """Print formatted hindcast report."""
    print(f"\n{title}")
    print("=" * 75)
    print(f"{'System':<14s} {'R²':>7s}  {'RMSE':>10s}  "
          f"{'NRMSE':>6s}  {'Bias':>10s}  {'Grade':>5s}")
    print("-" * 75)

    total_r2 = 0
    for name, score in scores.items():
        print(f"{name:<14s} {score.r_squared:>+7.3f}  "
              f"{score.rmse:>7.2f} {score.units:<3s}  "
              f"{score.nrmse:>6.2f}  "
              f"{score.bias:>+7.2f} {score.units:<3s}  "
              f"{score.grade():>5s}")
        total_r2 += score.r_squared

    if scores:
        mean_r2 = total_r2 / len(scores)
        print("-" * 75)
        print(f"{'Mean R²':<14s} {mean_r2:>+7.3f}")

    print()
    # Per-year detail for each system
    for name, score in scores.items():
        print(f"\n  {name} ({score.units}):")
        print(f"  {'Year':>6s}  {'Observed':>10s}  {'Predicted':>10s}  "
              f"{'Error':>10s}")
        for i, yr in enumerate(score.years):
            err = score.predicted[i] - score.observed[i]
            print(f"  {yr:>6d}  {score.observed[i]:>10.2f}  "
                  f"{score.predicted[i]:>10.2f}  {err:>+10.2f}")


def plot_hindcast(scores: Dict[str, HindcastScore],
                  save_path: str = None):
    """Plot observed vs predicted for all subsystems."""
    import matplotlib.pyplot as plt

    n = len(scores)
    fig, axes = plt.subplots(n, 1, figsize=(10, 3.5 * n), sharex=True)
    if n == 1:
        axes = [axes]

    for ax, (name, score) in zip(axes, scores.items()):
        ax.plot(score.years, score.observed, 'ko-', label='Observed',
                markersize=5)
        ax.plot(score.years, score.predicted, 'r^--', label='Model',
                markersize=5, alpha=0.8)
        ax.set_ylabel(f'{SCALES[name].observable}\n({score.units})',
                       fontsize=10)
        ax.legend(fontsize=9)
        ax.set_title(f'{name.capitalize()}: R²={score.r_squared:.3f}, '
                     f'RMSE={score.rmse:.2f} {score.units}, '
                     f'Grade={score.grade()}',
                     fontsize=11)
        ax.grid(alpha=0.3)

    axes[-1].set_xlabel('Year', fontsize=11)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='CEED Hindcast Validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python simulation/hindcast.py                     # Full report
    python simulation/hindcast.py --system solar      # Solar only
    python simulation/hindcast.py --compare           # Baseline vs anthro
    python simulation/hindcast.py --plot              # With plots
    python simulation/hindcast.py --start 2015        # Start from 2015
""")
    parser.add_argument('--system', type=str, default=None,
                        choices=SYSTEMS,
                        help='Run hindcast for a single subsystem')
    parser.add_argument('--compare', action='store_true',
                        help='Compare baseline vs anthropogenic')
    parser.add_argument('--start', type=int, default=2010,
                        help='Start year (default: 2010)')
    parser.add_argument('--end', type=int, default=2024,
                        help='End year (default: 2024)')
    parser.add_argument('--plot', action='store_true',
                        help='Show plots')
    parser.add_argument('--save', type=str, default=None,
                        help='Save plot to file')

    args = parser.parse_args()
    systems = [args.system] if args.system else None

    if args.compare:
        base_scores, anthro_scores = run_comparison(
            start_year=args.start, end_year=args.end)
        print_report(base_scores, "BASELINE (no anthropogenic forcing)")
        print_report(anthro_scores, "WITH ANTHROPOGENIC FORCING")

        print("\n" + "=" * 75)
        print("COMPARISON: R² improvement from anthropogenic forcing")
        print("-" * 75)
        for name in base_scores:
            if name in anthro_scores:
                delta = anthro_scores[name].r_squared - base_scores[name].r_squared
                print(f"  {name:<14s}: {delta:+.3f} "
                      f"({base_scores[name].grade()} -> {anthro_scores[name].grade()})")
    else:
        params = SystemParameters(anthropogenic=AnthropogenicForcing())
        scores = run_hindcast(params, systems=systems,
                              start_year=args.start, end_year=args.end)
        print_report(scores, f"Hindcast {args.start}-{args.end}")

        if args.plot or args.save:
            plot_hindcast(scores, save_path=args.save)


if __name__ == "__main__":
    main()
