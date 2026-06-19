#!/usr/bin/env python3
"""Compare gelation time from LigMC simulations against MD reference.

Uses trial results stored by main_ligmc.py. Outputs a CSV and a parity plot:
- Scatter of obtained (LigMC) gelation time vs expected (MD) gelation time.
- Diagonal y=x line highlights deviation from perfect agreement.
"""
from __future__ import annotations

import argparse
import os
import pickle
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class GelRecord:
    L: int
    mring: int
    nring: int
    mlin: int
    nlin: int
    md_t: float
    md_t_std: float | None = None


def _system_dir(root: str, rec: GelRecord) -> str:
    name = f"L{rec.L}_mring{rec.mring}_nring{rec.nring}_mlin{rec.mlin}_nlin{rec.nlin}"
    return os.path.join(root, name)


def _load_trial_results(run_dir: str) -> list[dict]:
    results_all = os.path.join(run_dir, "results_all.pkl")
    if os.path.exists(results_all):
        with open(results_all, "rb") as handle:
            payload = pickle.load(handle)
        if isinstance(payload, dict):
            for value in payload.values():
                return value
    for fname in os.listdir(run_dir):
        if fname.startswith("results_") and fname.endswith(".pkl"):
            with open(os.path.join(run_dir, fname), "rb") as handle:
                return pickle.load(handle)
    raise FileNotFoundError(f"No results_* or results_all.pkl found in {run_dir}")


def _extract_stages_to_half(trial_results: list[dict]) -> np.ndarray:
    values = []
    for item in trial_results:
        val = item.get("stages_to_half")
        if val is None:
            values.append(np.nan)
        else:
            values.append(float(val))
    return np.asarray(values, dtype=float)


def load_gel_point_file(path: str) -> list[GelRecord]:
    df = pd.read_csv(path)
    has_std = "t_std" in df.columns
    records: list[GelRecord] = []
    for _, row in df.iterrows():
        records.append(
            GelRecord(
                L=int(row["L"]),
                mring=int(row["Mring"]),
                nring=int(row["Nring"]),
                mlin=int(row["Mlin"]),
                nlin=int(row["Nlin"]),
                md_t=float(row["t"]),
                md_t_std=float(row["t_std"]) if has_std else None,
            )
        )
    return records


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--gel-file", required=True, help="Path to gel_point_time.txt")
    p.add_argument("--results-root", required=True, help="Root directory containing per-system output folders")
    p.add_argument("--out-csv", default="gel_time_compare.csv")
    p.add_argument("--out-png", default="gel_time_compare.png")
    args = p.parse_args()

    records = load_gel_point_file(args.gel_file)

    rows = []
    missing = []
    for rec in records:
        run_dir = _system_dir(args.results_root, rec)
        if not os.path.isdir(run_dir):
            missing.append(run_dir)
            continue
        try:
            trials = _load_trial_results(run_dir)
        except FileNotFoundError:
            missing.append(run_dir)
            continue

        stages = _extract_stages_to_half(trials)
        mean = float(np.nanmean(stages)) if np.any(np.isfinite(stages)) else float("nan")
        std = float(np.nanstd(stages, ddof=1)) if np.sum(np.isfinite(stages)) > 1 else 0.0
        rows.append(
            {
                "L": rec.L,
                "Mring": rec.mring,
                "Nring": rec.nring,
                "Mlin": rec.mlin,
                "Nlin": rec.nlin,
                "md_t": rec.md_t,
                "md_t_std": float(rec.md_t_std) if rec.md_t_std is not None else np.nan,
                "sim_mean": mean,
                "sim_std": std,
                "n_trials": int(np.sum(np.isfinite(stages))),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(args.out_csv, index=False)

    if df.empty:
        raise SystemExit("No matching results found. Check results-root and run directories.")

    # --- extract arrays ------------------------------------------------
    expected = df["md_t"].to_numpy(dtype=float)            # MD reference
    expected_err = df["md_t_std"].to_numpy(dtype=float) if "md_t_std" in df.columns else None
    obtained = df["sim_mean"].to_numpy(dtype=float)        # LigMC result
    obtained_err = df["sim_std"].to_numpy(dtype=float)

    finite = np.isfinite(expected) & np.isfinite(obtained)
    expected = expected[finite]
    obtained = obtained[finite]
    obtained_err = obtained_err[finite]
    if expected_err is not None:
        expected_err = expected_err[finite]

    mae = float(np.mean(np.abs(obtained - expected)))
    rmse = float(np.sqrt(np.mean((obtained - expected) ** 2)))

    # --- parity plot ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 6))

    ax.errorbar(
        expected,
        obtained,
        xerr=expected_err,
        yerr=obtained_err,
        fmt="o",
        color="#2563eb",
        markeredgecolor="white",
        markeredgewidth=0.6,
        markersize=7,
        ecolor="#94a3b8",
        elinewidth=1.2,
        capsize=3,
        alpha=0.9,
        zorder=3,
        label="Systems",
    )

    # diagonal y = x line with padding
    lo = float(min(np.min(expected), np.min(obtained)))
    hi = float(max(np.max(expected), np.max(obtained)))
    pad = 0.10 * (hi - lo) if hi > lo else 1.0
    diag = [lo - pad, hi + pad]
    ax.plot(diag, diag, ls="--", color="#dc2626", lw=1.4, alpha=0.7, zorder=2,
            label="Perfect agreement")

    ax.set_xlim(diag)
    ax.set_ylim(diag)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("Expected gelation time — MD (stages)", fontsize=11)
    ax.set_ylabel("Obtained gelation time — LigMC (stages)", fontsize=11)
    ax.set_title(f"Gel-point parity  (MAE = {mae:.2f},  RMSE = {rmse:.2f})",
                 fontsize=12, fontweight="semibold")
    ax.legend(fontsize=9, framealpha=0.8)
    ax.grid(alpha=0.20, linewidth=0.6)

    fig.tight_layout()
    fig.savefig(args.out_png, dpi=160)

    if missing:
        print("Missing results for:")
        for path in missing:
            print(f"- {path}")

    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.out_png}")


if __name__ == "__main__":
    main()
