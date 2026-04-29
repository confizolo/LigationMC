"""Parameter sweep script for LigMC gelation across nring and nlin grids."""

from __future__ import annotations

import argparse
import multiprocessing
import os
import sys
import traceback
from typing import Any

import numpy as np

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None

from .analysis import save_json, save_results_all
from .main_ligmc import FITTED_K1_DEFAULT, FITTED_K2_DEFAULT, FITTED_A_DEFAULT, _first_stage_to_half, run_single_trial
from .polymer_utils import RESULTS_DIR, calculate_polymer_numbers


def build_system_grid(
    nrings: np.ndarray,
    nlins: np.ndarray,
    ccsr: float = 5.0,
    ccsl: float = 0.05,
    L: float = 80.0,
) -> list[dict[str, Any]]:
    """Build list of system config dicts from the (nring, nlin) grid.

    For each (nring, nlin) pair, computes (mring, mlin) via calculate_polymer_numbers.
    """
    systems = []
    for nring in nrings:
        for nlin in nlins:
            mring, mlin = calculate_polymer_numbers(ccsr, ccsl, Nr=int(nring), Nl=int(nlin), L=L)
            systems.append(
                {
                    "L": L,
                    "mring": mring,
                    "nring": int(nring),
                    "mlin": mlin,
                    "nlin": int(nlin),
                }
            )
    return systems


def _run_single_trial_safe(config: dict[str, Any]) -> dict[str, Any] | None:
    try:
        return run_single_trial(config)
    except Exception as e:
        print(f"\n[ERROR] Trial failed for seed={config.get('seed')}: {e}")
        traceback.print_exc(file=sys.stdout)
        return None


def run_trials_for_system(
    sys_cfg: dict[str, Any],
    trials: int,
    n_stages: int,
    workers: int,
    k1: float,
    k2: float,
    alpha: float,
    nu: float,
    val_A: float,
    max_steps: int,
) -> list[dict[str, Any]]:
    """Run a specified number of trials for a single system configuration."""
    configs = []
    for t in range(trials):
        seed = hash((sys_cfg["nring"], sys_cfg["nlin"], t)) % (2**32)
        configs.append(
            {
                **sys_cfg,
                "n_stages": n_stages,
                "k1": k1,
                "k2": k2,
                "alpha": alpha,
                "nu": nu,
                "val_A": val_A,
                "max_steps": max_steps,
                "seed": seed,
            }
        )

    results = []
    if workers > 1:
        with multiprocessing.Pool(processes=workers) as pool:
            for res in pool.imap_unordered(_run_single_trial_safe, configs):
                if res is not None:
                    results.append(res)
    else:
        for cfg in configs:
            res = _run_single_trial_safe(cfg)
            if res is not None:
                results.append(res)

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Sweep LigMC gelation over nring/nlin grid.")
    parser.add_argument("--nring_min", type=int, default=256)
    parser.add_argument("--nring_max", type=int, default=1040)
    parser.add_argument("--nring_step", type=int, default=16)
    parser.add_argument("--nlin_min", type=int, default=64)
    parser.add_argument("--nlin_max", type=int, default=224)
    parser.add_argument("--nlin_step", type=int, default=16)
    parser.add_argument("--ccsr", type=float, default=5.0)
    parser.add_argument("--ccsl", type=float, default=0.05)
    parser.add_argument("--L", type=float, default=80.0)
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--n_stages", type=int, default=100)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--skip_plots", action="store_true", help="Skip plotting step")
    parser.add_argument("--out_dir", type=str, default=RESULTS_DIR)
    
    # Engine defaults matching main_ligmc
    parser.add_argument("--k1", type=float, default=FITTED_K1_DEFAULT)
    parser.add_argument("--k2", type=float, default=FITTED_K2_DEFAULT)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--nu", type=float, default=0.5)
    parser.add_argument("--val_A", type=float, default=FITTED_A_DEFAULT)
    parser.add_argument("--max_steps", type=int, default=1000000)

    args = parser.parse_args()

    nrings = np.arange(args.nring_min, args.nring_max, args.nring_step)
    nlins = np.arange(args.nlin_min, args.nlin_max, args.nlin_step)
    
    systems = build_system_grid(nrings, nlins, args.ccsr, args.ccsl, args.L)
    print(f"Built grid of {len(systems)} systems.")

    os.makedirs(args.out_dir, exist_ok=True)
    csv_path = os.path.join(args.out_dir, "sweep_summary.csv")
    
    # Write header for the summary CSV if it doesn't exist
    if not os.path.exists(csv_path):
        with open(csv_path, "w") as f:
            f.write("L,mring,nring,mlin,nlin,mean_stages,std_stages,n_trials\n")

    sweep_results = {}
    
    iterator = tqdm(systems, desc="Systems") if tqdm else systems
    
    for sys_cfg in iterator:
        tag = f"L{sys_cfg['L']}_mring{sys_cfg['mring']}_nring{sys_cfg['nring']}_mlin{sys_cfg['mlin']}_nlin{sys_cfg['nlin']}"
        
        trial_results = run_trials_for_system(
            sys_cfg=sys_cfg,
            trials=args.trials,
            n_stages=args.n_stages,
            workers=args.workers,
            k1=args.k1,
            k2=args.k2,
            alpha=args.alpha,
            nu=args.nu,
            val_A=args.val_A,
            max_steps=args.max_steps,
        )
        
        if not trial_results:
            print(f"WARNING: All trials failed for {tag}")
            continue
            
        stages_to_half_values = [t["stages_to_half"] for t in trial_results if t["stages_to_half"] is not None]
        
        if stages_to_half_values:
            mean_s = float(np.mean(stages_to_half_values))
            std_s = float(np.std(stages_to_half_values))
        else:
            mean_s = np.nan
            std_s = np.nan
            
        sys_result = {
            **sys_cfg,
            "mean_stages": mean_s,
            "std_stages": std_s,
            "n_trials_successful": len(trial_results),
            "n_gelled": len(stages_to_half_values)
        }
        
        sweep_results[tag] = sys_result
        
        # Append to CSV
        with open(csv_path, "a") as f:
            f.write(f"{sys_cfg['L']},{sys_cfg['mring']},{sys_cfg['nring']},{sys_cfg['mlin']},{sys_cfg['nlin']},"
                    f"{mean_s:.4f},{std_s:.4f},{len(stages_to_half_values)}\n")

    save_results_all(sweep_results, out_dir=args.out_dir)
    print(f"\nSweep complete. Summary saved to {csv_path}")

if __name__ == "__main__":
    main()
