"""Analysis and plotting helpers for LigMC outputs."""

from __future__ import annotations

import json
import os
import pickle
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from polymer_utils import RESULTS_DIR, calculate_valence


def _system_tag(L: int, mring: int, nring: int, mlin: int, nlin: int) -> str:
    return f"L{L}_mring{mring}_nring{nring}_mlin{mlin}_nlin{nlin}"


def save_trial_results(
    trial_results: list[dict[str, Any]],
    L: int,
    mring: int,
    nring: int,
    mlin: int,
    nlin: int,
    out_dir: str = RESULTS_DIR,
) -> str:
    os.makedirs(out_dir, exist_ok=True)
    filename = f"results_{_system_tag(L, mring, nring, mlin, nlin)}.pkl"
    path = os.path.join(out_dir, filename)
    with open(path, "wb") as handle:
        pickle.dump(trial_results, handle)
    return path


def save_results_all(results_by_system: dict[str, Any], out_dir: str = RESULTS_DIR) -> str:
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "results_all.pkl")
    with open(path, "wb") as handle:
        pickle.dump(results_by_system, handle)
    return path


def save_event_timeline(event_timeline: list[dict[str, Any]], out_dir: str = RESULTS_DIR) -> str:
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "event_timeline.pkl")
    with open(path, "wb") as handle:
        pickle.dump(event_timeline, handle)
    return path


def save_json(payload: dict[str, Any], filename: str, out_dir: str = RESULTS_DIR) -> str:
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, filename)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    return path


def plot_gelation_curves(trial_results: list[dict[str, Any]], out_dir: str = RESULTS_DIR) -> str:
    os.makedirs(out_dir, exist_ok=True)

    if not trial_results:
        raise ValueError("No trial results available for gelation-curve plotting.")

    max_len = max(len(item["largest_component_fraction"]) for item in trial_results)
    matrix = np.full((len(trial_results), max_len), np.nan, dtype=float)

    for i, item in enumerate(trial_results):
        arr = np.asarray(item["largest_component_fraction"], dtype=float)
        matrix[i, : len(arr)] = arr

    mean_curve = np.nanmean(matrix, axis=0)
    low = np.nanpercentile(matrix, 10, axis=0)
    high = np.nanpercentile(matrix, 90, axis=0)

    x = np.arange(1, max_len + 1)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, mean_curve, lw=2, label="mean")
    ax.fill_between(x, low, high, alpha=0.25, label="10-90 percentile")
    ax.axhline(0.5, color="k", ls="--", lw=1)
    ax.set_xlabel("Stage")
    ax.set_ylabel("Largest component fraction")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()

    path = os.path.join(out_dir, "gelation_curves.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_degree_distributions(
    degree_distributions: list[dict[int, int]],
    out_dir: str = RESULTS_DIR,
) -> str:
    os.makedirs(out_dir, exist_ok=True)
    if not degree_distributions:
        raise ValueError("No degree distributions provided.")

    n = len(degree_distributions)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)

    for idx, dist in enumerate(degree_distributions):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        keys = sorted(dist)
        vals = np.array([dist[k] for k in keys], dtype=float)
        vals = vals / vals.sum() if vals.sum() > 0 else vals
        ax.bar(keys, vals, alpha=0.8, edgecolor="k")
        ax.set_xlabel("Degree")
        ax.set_ylabel("Fraction")
        ax.set_title(f"Trial {idx + 1}")
        ax.grid(alpha=0.2)

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].axis("off")

    fig.tight_layout()
    path = os.path.join(out_dir, "degree_distributions.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_stages_vs_nlin_nring(records: list[dict[str, float]], out_dir: str = RESULTS_DIR) -> str:
    os.makedirs(out_dir, exist_ok=True)
    if not records:
        raise ValueError("No records provided for stage-parameter plot.")

    nlin = np.array([r["nlin"] for r in records], dtype=float)
    nring = np.array([r["nring"] for r in records], dtype=float)
    mlin = np.array([r["mlin"] for r in records], dtype=float)
    stages_raw = np.array([r["stages_to_half"] for r in records], dtype=float)
    valid = np.isfinite(stages_raw) & (stages_raw > 0)
    if not np.any(valid):
        raise ValueError("No positive finite stages_to_half values available for log-scale plotting.")

    nlin = nlin[valid]
    nring = nring[valid]
    mlin = mlin[valid]
    stages = stages_raw[valid]

    valence = np.array([calculate_valence(int(nl), int(nr)) for nl, nr in zip(nlin, nring)])

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    sc = axes[0, 0].scatter(valence * mlin, stages, c=nlin, cmap="viridis", edgecolor="k")
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel("valence * mlin")
    axes[0, 0].set_ylabel("stages to 50%")
    fig.colorbar(sc, ax=axes[0, 0], label="nlin")

    axes[0, 1].scatter(nlin, stages, c=nring, cmap="plasma", edgecolor="k")
    axes[0, 1].set_xlabel("nlin")
    axes[0, 1].set_ylabel("stages to 50%")
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_yscale("log")

    axes[1, 0].scatter(nring, stages, c=nlin, cmap="magma", edgecolor="k")
    axes[1, 0].set_xlabel("nring")
    axes[1, 0].set_ylabel("stages to 50%")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_yscale("log")

    axes[1, 1].scatter(nring / nlin, stages, c=valence, cmap="cividis", edgecolor="k")
    axes[1, 1].set_xlabel("nring/nlin")
    axes[1, 1].set_ylabel("stages to 50%")
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_yscale("log")

    for ax in axes.flat:
        ax.grid(alpha=0.2)

    fig.tight_layout()
    path = os.path.join(out_dir, "stages_vs_nlin_nring.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
