"""Main entry point for staged LigMC simulations."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from dataclasses import asdict
import sys
from typing import Any

# Ensure repo root is on sys.path for direct script execution.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import numpy as np

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None

from simulation.analysis import (
    plot_degree_distributions,
    plot_gelation_curves,
    plot_stages_vs_nlin_nring,
    save_event_timeline,
    save_results_all,
    save_trial_results,
)
from simulation.dsmc_engine import CyclisationEvent, DSMCEngine
from simulation.network_builder import NetworkBuilder
from simulation.polymer_utils import RESULTS_DIR


FITTED_K1_DEFAULT = 1.0
FITTED_K2_DEFAULT = 12840.849325710129
FITTED_A_DEFAULT = 0.20927677484111143


def _first_stage_to_half(fracs: list[float]) -> int | None:
    for idx, val in enumerate(fracs, start=1):
        if val >= 0.5:
            return idx
    return None


def _run_stage(
    network: NetworkBuilder,
    config: dict[str, Any],
    stage: int,
    phi_ref: float,
    rng: np.random.Generator,
) -> tuple[float, float, int, list[dict[str, Any]]]:
    """Run one stage: inject linears, DSMC until exhausted, process events.

    Returns (largest_component_frac, box_length, n_events, event_records).
    """
    total_ring_monomers = float(sum(network.ring_lengths.values()))
    total_monomers = total_ring_monomers + float(config["mlin"] * config["nlin"])
    box_volume = total_monomers / phi_ref
    stage_L = float(box_volume ** (1.0 / 3.0))

    linears = [config["nlin"]] * config["mlin"]
    engine = DSMCEngine(
        linear_lengths=linears,
        k1=config["k1"],
        k2=config["k2"],
        alpha=config["alpha"],
        nu=config["nu"],
        rng=rng,
    )

    events = engine.run_until_exhausted(max_steps=config["max_steps"])
    n_events = len(events)
    event_records = []

    for event in events:
        if isinstance(event, CyclisationEvent):
            event = network.process_cyclisation(
                event,
                A=config["val_A"],
                box_volume=box_volume,
            )
            event_records.append(
                {
                    "stage": stage,
                    "L": stage_L,
                    "event_type": "cyclisation",
                    **asdict(event),
                }
            )
        else:
            event_records.append(
                {
                    "stage": stage,
                    "L": stage_L,
                    "event_type": "merge",
                    **asdict(event),
                }
            )

    frac = network.largest_component_fraction()
    return frac, stage_L, n_events, event_records


def run_single_trial(config: dict[str, Any]) -> dict[str, Any]:
    rng = np.random.default_rng(config["seed"])

    network = NetworkBuilder(
        initial_ring_lengths=[config["nring"]] * config["mring"],
        rng=rng,
    )

    stage_fractions: list[float] = []
    stage_events: list[int] = []
    stage_box_lengths: list[float] = []
    degree_distribution: dict[int, int] = {}
    event_timeline: list[dict[str, Any]] = []

    # Keep volume fraction constant across stages by rescaling L with total monomers.
    # Reference volume fraction is based on initial rings + first-stage linears.
    initial_monomers = float(config["mring"] * config["nring"] + config["mlin"] * config["nlin"])
    if initial_monomers <= 0.0:
        raise ValueError("Total monomers must be positive to set volume fraction.")
    phi_ref = initial_monomers / float(config["L"] ** 3)

    for stage in range(config["n_stages"]):
        frac, stage_L, n_events, event_records = _run_stage(network, config, stage, phi_ref, rng)
        stage_fractions.append(frac)
        stage_box_lengths.append(stage_L)
        stage_events.append(n_events)
        event_timeline.extend(event_records)
        degree_distribution = network.degree_distribution()

    return {
        "seed": config["seed"],
        "largest_component_fraction": stage_fractions,
        "stages_to_half": _first_stage_to_half(stage_fractions),
        "degree_distribution": degree_distribution,
        "events_per_stage": stage_events,
        "box_length_per_stage": stage_box_lengths,
        "final_n_rings": int(network.graph.number_of_nodes()),
        "event_timeline": event_timeline,
    }


def _run_trials(args: argparse.Namespace) -> list[dict[str, Any]]:
    base_seed = args.seed if args.seed is not None else int(np.random.SeedSequence().entropy)
    seed_seq = np.random.SeedSequence(base_seed)
    trial_seeds = [int(s.generate_state(1)[0]) for s in seed_seq.spawn(args.trials)]

    jobs = []
    for trial_idx, trial_seed in enumerate(trial_seeds):
        jobs.append(
            {
                "trial_idx": trial_idx,
                "seed": trial_seed,
                "L": args.L,
                "mring": args.mring,
                "nring": args.nring,
                "mlin": args.mlin,
                "nlin": args.nlin,
                "n_stages": args.n_stages,
                "k1": args.k1,
                "k2": args.k2,
                "alpha": args.alpha,
                "nu": args.nu,
                "val_A": args.val_A,
                "max_steps": args.max_steps,
            }
        )

    def _progress(iterable, enabled: bool, total: int) -> Any:
        if enabled and tqdm is not None:
            return tqdm(iterable, total=total, desc="Trials")
        return iterable

    if args.workers > 1:
        with mp.Pool(processes=args.workers) as pool:
            if args.progress and tqdm is not None:
                results = list(_progress(pool.imap(run_single_trial, jobs), True, total=len(jobs)))
            else:
                results = pool.map(run_single_trial, jobs)
    else:
        results = [run_single_trial(job) for job in _progress(jobs, args.progress, total=len(jobs))]

    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="LigMC staged DSMC/network-growth simulation")
    parser.add_argument("--L", type=int, default=80, help="Simulation box length label for output naming")
    parser.add_argument("--mring", type=int, default=222, help="Initial number of ring polymers")
    parser.add_argument("--nring", type=int, default=256, help="Ring polymer length")
    parser.add_argument("--mlin", type=int, default=2, help="Linears injected per stage")
    parser.add_argument("--nlin", type=int, default=128, help="Linear polymer length")
    parser.add_argument("--n_stages", type=int, default=100, help="Number of stages")
    parser.add_argument("--trials", type=int, default=1000, help="Number of independent trials")
    parser.add_argument("--k1", type=float, default=FITTED_K1_DEFAULT, help="Merging rate constant")
    parser.add_argument("--k2", type=float, default=FITTED_K2_DEFAULT, help="Cyclisation rate constant")
    parser.add_argument("--alpha", type=float, default=1.0, help="Kernel exponent alpha")
    parser.add_argument("--nu", type=float, default=0.5, help="Kernel exponent nu")
    parser.add_argument("--val_A", type=float, default=FITTED_A_DEFAULT, help="Valence model prefactor A")
    parser.add_argument("--max_steps", type=int, default=50000, help="Max SSA events per stage")
    parser.add_argument("--seed", type=int, default=None, help="Base RNG seed")
    parser.add_argument("--workers", type=int, default=1, help="Multiprocessing workers")
    parser.add_argument("--progress", action="store_true", help="Show a progress bar (requires tqdm)")
    parser.add_argument("--out_dir", type=str, default=RESULTS_DIR, help="Output directory")
    parser.add_argument("--skip_plots", action="store_true", help="Skip plotting outputs")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    system_dir_name = f"L{args.L}_mring{args.mring}_nring{args.nring}_mlin{args.mlin}_nlin{args.nlin}"
    run_out_dir = os.path.join(args.out_dir, system_dir_name)
    os.makedirs(run_out_dir, exist_ok=True)
    trial_results = _run_trials(args)

    save_trial_results(
        trial_results,
        L=args.L,
        mring=args.mring,
        nring=args.nring,
        mlin=args.mlin,
        nlin=args.nlin,
        out_dir=run_out_dir,
    )

    system_key = f"{args.L} {args.mring} {args.nring} {args.mlin} {args.nlin}"
    save_results_all({system_key: trial_results}, out_dir=run_out_dir)

    event_timeline: list[dict[str, Any]] = []
    for item in trial_results:
        event_timeline.extend(item["event_timeline"])
    save_event_timeline(event_timeline, out_dir=run_out_dir)

    if not args.skip_plots:
        plot_gelation_curves(trial_results, out_dir=run_out_dir)

        stages_records = [
            {
                "nlin": float(args.nlin),
                "nring": float(args.nring),
                "mlin": float(args.mlin),
                "stages_to_half": float(item["stages_to_half"] or np.nan),
            }
            for item in trial_results
        ]
        try:
            plot_stages_vs_nlin_nring(stages_records, out_dir=run_out_dir)
        except ValueError:
            # Small/sparse runs may never reach 50% and cannot be log-plotted.
            pass


if __name__ == "__main__":
    main()
