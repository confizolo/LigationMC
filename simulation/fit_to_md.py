"""Fit k1/k2 from cyclized-length PMF using Gillespie SSA."""

from __future__ import annotations

import argparse
import math
import os
from collections import Counter
import sys
from typing import Any

# Ensure repo root is on sys.path for direct script execution.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import matplotlib.pyplot as plt
import numpy as np

from simulation.analysis import save_json
from simulation.dsmc_engine import CyclisationEvent, DSMCEngine
from simulation.polymer_utils import RESULTS_DIR

try:
    from scipy.optimize import minimize
except ImportError as exc:  # pragma: no cover
    raise ImportError("scipy is required for fit_to_md.py: pip install scipy") from exc


TARGET_PMF: dict[int, float] = {
    64: 0.7799,
    128: 0.1557,
    192: 0.0371,
    256: 0.0222,
    320: 0.0052,
    384: 0.0015,
}


def _normalize_pmf(pmf: dict[int, float]) -> dict[int, float]:
    total = float(sum(max(v, 0.0) for v in pmf.values()))
    if total <= 0.0:
        return {}
    return {int(k): float(max(v, 0.0) / total) for k, v in pmf.items() if v > 0.0}


def _js_divergence(p: dict[int, float], q: dict[int, float]) -> float:
    eps = 1e-12
    keys = sorted(set(p) | set(q))
    p_arr = np.array([p.get(k, 0.0) for k in keys], dtype=float)
    q_arr = np.array([q.get(k, 0.0) for k in keys], dtype=float)
    p_arr = np.clip(p_arr, eps, None)
    q_arr = np.clip(q_arr, eps, None)
    p_arr /= float(p_arr.sum())
    q_arr /= float(q_arr.sum())
    m_arr = 0.5 * (p_arr + q_arr)

    kl_pm = float(np.sum(p_arr * np.log(p_arr / m_arr)))
    kl_qm = float(np.sum(q_arr * np.log(q_arr / m_arr)))
    return 0.5 * (kl_pm + kl_qm)


def simulate_length_pmf(
    k1: float,
    k2: float,
    nlin: int,
    mlin: int,
    nu: float,
    n_trials: int,
    rng: np.random.Generator,
    max_steps: int = 50000,
) -> dict[int, float]:
    """Run n_trials SSA stages and return normalized cyclized-length PMF."""
    counts: Counter[int] = Counter()

    for _ in range(n_trials):
        linears = [nlin] * mlin
        engine = DSMCEngine(linear_lengths=linears, k1=k1, k2=k2, nu=nu, rng=rng)
        events = engine.run_until_exhausted(max_steps=max_steps)
        for event in events:
            if isinstance(event, CyclisationEvent):
                counts[int(event.ring_length)] += 1

    total = int(sum(counts.values()))
    if total <= 0:
        return {}
    return {length: count / total for length, count in sorted(counts.items())}


def fit_k1_k2(
    target_pmf: dict[int, float],
    nlin: int = 64,
    mlin: int = 6,
    nu: float = 0.5,
    n_trials: int = 5000,
    k1_fixed: float = 1.0,
    fit_both: bool = False,
    seed: int | None = None,
) -> dict[str, float]:
    """Fit k1, k2 by minimizing JS-divergence to target PMF.

    Returns {'k1': float, 'k2': float, 'kappa': float, 'js_div': float}.
    """
    if seed is None:
        seed = 1234
    target = _normalize_pmf(target_pmf)

    if fit_both:

        def objective(x: np.ndarray) -> float:
            k1 = float(10.0 ** x[0])
            k2 = float(10.0 ** x[1])
            sim_pmf = simulate_length_pmf(k1, k2, nlin, mlin, nu, n_trials, rng)
            return _js_divergence(sim_pmf, target)

        x0 = np.array([math.log10(max(k1_fixed, 1e-12)), -2.0], dtype=float)
        opt = minimize(objective, x0=x0, method="Nelder-Mead")
        k1 = float(10.0 ** opt.x[0])
        k2 = float(10.0 ** opt.x[1])
        js_div = float(opt.fun)
    else:

        cache: dict[float, float] = {}

        def _objective_from_log10_k2(log10_k2: float) -> float:
            key = float(round(log10_k2, 6))
            if key in cache:
                return cache[key]

            k2 = float(10.0 ** log10_k2)
            # Re-seed every evaluation so the optimizer sees a stable objective.
            rng_eval = np.random.default_rng(seed)
            sim_pmf = simulate_length_pmf(k1_fixed, k2, nlin, mlin, nu, n_trials, rng_eval)
            value = _js_divergence(sim_pmf, target)
            cache[key] = value
            return value

        def objective_1d(x: np.ndarray) -> float:
            return _objective_from_log10_k2(float(x[0]))

        # Coarse scan in log-space improves robustness before local simplex refinement.
        grid = np.linspace(-8.0, 4.0, 25)
        grid_vals = np.array([_objective_from_log10_k2(float(g)) for g in grid], dtype=float)
        x0 = np.array([float(grid[int(np.argmin(grid_vals))])], dtype=float)
        opt = minimize(objective_1d, x0=x0, method="Nelder-Mead")
        k1 = float(k1_fixed)
        k2 = float(10.0 ** opt.x[0])
        js_div = float(opt.fun)

    kappa = float(k2 / k1) if k1 > 0 else float("inf")
    return {"k1": k1, "k2": k2, "kappa": kappa, "js_div": js_div}


def _plot_pmf_comparison(
    target_pmf: dict[int, float],
    simulated_pmf: dict[int, float],
    out_dir: str,
) -> str:
    keys = sorted(set(target_pmf) | set(simulated_pmf))
    x = np.arange(len(keys), dtype=float)
    width = 0.42

    y_target = np.array([target_pmf.get(k, 0.0) for k in keys], dtype=float)
    y_sim = np.array([simulated_pmf.get(k, 0.0) for k in keys], dtype=float)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(x - width / 2.0, y_target, width=width, label="MD target", alpha=0.85)
    ax.bar(x + width / 2.0, y_sim, width=width, label="SSA fit", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([str(k) for k in keys])
    ax.set_xlabel("Cyclized length")
    ax.set_ylabel("PMF")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "fit_to_md_pmf_comparison.png")
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit k1/k2 by matching the cyclized-length PMF.")
    parser.add_argument("--nlin", type=int, default=64)
    parser.add_argument("--mlin", type=int, default=6)
    parser.add_argument("--nu", type=float, default=0.5)
    parser.add_argument("--n_trials", type=int, default=5000)
    parser.add_argument("--k1_fixed", type=float, default=1.0, help="Reference k1 when fitting only k2")
    parser.add_argument("--fit_both", action="store_true", help="Fit k1 and k2 jointly (2D search)")
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--out_dir", type=str, default=RESULTS_DIR)
    args = parser.parse_args()

    fit = fit_k1_k2(
        target_pmf=TARGET_PMF,
        nlin=args.nlin,
        mlin=args.mlin,
        nu=args.nu,
        n_trials=args.n_trials,
        k1_fixed=args.k1_fixed,
        fit_both=args.fit_both,
        seed=args.seed,
    )

    rng = np.random.default_rng(args.seed + 1)
    sim_pmf = simulate_length_pmf(
        fit["k1"],
        fit["k2"],
        nlin=args.nlin,
        mlin=args.mlin,
        nu=args.nu,
        n_trials=args.n_trials,
        rng=rng,
    )

    fit["target_pmf"] = {str(k): float(v) for k, v in TARGET_PMF.items()}
    fit["simulated_pmf"] = {str(k): float(v) for k, v in sim_pmf.items()}

    os.makedirs(args.out_dir, exist_ok=True)
    out_json = save_json(fit, filename="fitted_k1_k2.json", out_dir=args.out_dir)
    out_plot = _plot_pmf_comparison(TARGET_PMF, sim_pmf, out_dir=args.out_dir)

    print(f"k1      = {fit['k1']:.6e}")
    print(f"k2      = {fit['k2']:.6e}")
    print(f"kappa   = {fit['kappa']:.6e}")
    print(f"JS div  = {fit['js_div']:.6e}")
    print(f"saved   = {out_json}")
    print(f"plot    = {out_plot}")


if __name__ == "__main__":
    main()
