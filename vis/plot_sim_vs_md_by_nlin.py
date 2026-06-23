"""Create DSMC-vs-MD cyclized-length PMF overlays for multiple nlin values."""

from __future__ import annotations

import argparse
import json
import math
import os
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RESULTS_DIR = (
    "/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc"
)
FITTED_K1_DEFAULT = 1.0
FITTED_K2_DEFAULT = 7928.458898884091


def save_json(data: dict, filename: str, out_dir: str) -> str:
    """Write *data* as JSON to *out_dir*/*filename* and return the path."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, filename)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    return path


def _kl_divergence(p: dict[int, float], q: dict[int, float]) -> float:
    """KL(P || Q) over shared keys with p_i > 0 and q_i > 0."""
    kl = 0.0
    for k in p:
        pk = p[k]
        qk = q.get(k, 0.0)
        if pk > 0.0 and qk > 0.0:
            kl += pk * math.log(pk / qk)
    return kl


def _js_divergence(p: dict[int, float], q: dict[int, float]) -> float:
    """Jensen-Shannon divergence between two PMF dicts."""
    all_keys = set(p) | set(q)
    m: dict[int, float] = {}
    for k in all_keys:
        m[k] = 0.5 * (p.get(k, 0.0) + q.get(k, 0.0))
    return 0.5 * _kl_divergence(p, m) + 0.5 * _kl_divergence(q, m)


def simulate_dsmc_length_pmf(
    k1: float,
    k2: float,
    nlin: int,
    mlin: int,
    nu: float,
    n_trials: int,
    rng: np.random.Generator,
) -> dict[int, float]:
    """Particle DSMC producing the cyclised-ring-length PMF.

    Starts with *mlin* linears of length *nlin*.  Uses acceptance-rejection
    with majorant rates for merge (Smoluchowski kernel) and cyclisation.

    Returns a normalised PMF dict {length: probability}.
    """
    alpha = 1.0
    density = k1 * mlin
    counts: dict[int, int] = {}

    for _ in range(n_trials):
        masses = [nlin] * mlin
        n_chains = mlin

        # Initial majorant rate estimates (must be > 0)
        k_max = (nlin ** (-alpha) + nlin ** (-alpha)) * (nlin ** nu + nlin ** nu)
        r_max = k2 * (nlin ** (-4.0 * nu))
        if k_max <= 0:
            k_max = 1e-4
        if r_max <= 0:
            r_max = 1e-4

        while n_chains > 0:
            if n_chains > 1:
                p_ann = 1.0 / (
                    1.0
                    + (2.0 * mlin * r_max)
                    / ((n_chains - 1) * density * k_max)
                )
            else:
                p_ann = 0.0

            if rng.random() < p_ann:
                # Merge attempt
                active = [idx for idx, m in enumerate(masses) if m > 0]
                i = rng.choice(active)
                j = rng.integers(0, mlin)
                while masses[j] == 0 or j == i:
                    j = rng.integers(0, mlin)

                mi = masses[i]
                mj = masses[j]
                k_ij = (mi ** (-alpha) + mj ** (-alpha)) * (mi ** nu + mj ** nu)

                if k_ij > k_max:
                    k_max = k_ij
                else:
                    if rng.random() < k_ij / k_max:
                        masses[j] = mi + mj
                        masses[i] = 0
                        n_chains -= 1
            else:
                # Cyclisation attempt
                active = [idx for idx, m in enumerate(masses) if m > 0]
                k = rng.choice(active)
                mk = masses[k]
                rmk = k2 * (mk ** (-4.0 * nu))

                if rmk > r_max:
                    r_max = rmk
                else:
                    if rng.random() < rmk / r_max:
                        masses[k] = 0
                        n_chains -= 1
                        counts[mk] = counts.get(mk, 0) + 1

    total = sum(counts.values())
    if total <= 0:
        return {}
    return {k: v / total for k, v in sorted(counts.items())}


DEFAULT_MD_COUNTS_CSV = (
    "/storage/cmstore02/groups/TAPLab/fconforto-projects/"
    "fconforto-olympic-gels/results/histories/dist_cyclized_linear_length_by_nlin_all_systems.csv"
)
DEFAULT_FIT_JSON = os.path.join(RESULTS_DIR, "fitted_k1_k2.json")


def _parse_nlins(raw: str) -> list[int]:
    vals = [int(x.strip()) for x in raw.split(",") if x.strip()]
    uniq = sorted(set(vals))
    if not uniq:
        raise ValueError("No valid nlin values parsed from --nlins")
    return uniq


def _normalize_counts(sub: pd.DataFrame) -> dict[int, float]:
    if sub.empty:
        return {}
    counts = sub.groupby("linear_length", as_index=True)["count"].sum().to_dict()
    total = float(sum(max(v, 0.0) for v in counts.values()))
    if total <= 0:
        return {}
    return {int(k): float(v / total) for k, v in counts.items() if v > 0}


def _plot_overlay(
    md_pmf: dict[int, float],
    dsmc_pmf: dict[int, float],
    nlin: int,
    out_dir: str,
) -> str:
    keys = sorted(set(md_pmf) | set(dsmc_pmf))
    x = np.arange(len(keys), dtype=float)
    width = 0.38

    y_md = np.array([md_pmf.get(k, 0.0) for k in keys], dtype=float)
    y_dsmc = np.array([dsmc_pmf.get(k, 0.0) for k in keys], dtype=float)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(x - width / 2.0, y_md, width=width, label="MD target", alpha=0.85)
    ax.bar(x + width / 2.0, y_dsmc, width=width, label="DSMC", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([str(k) for k in keys])
    ax.set_xlabel("Cyclized length")
    ax.set_ylabel("PMF")
    ax.set_title(f"Cyclized-Length PMF Comparison (nlin={nlin})")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"fit_to_md_pmf_comparison_nlin{nlin}.png")
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate PMF overlays (DSMC vs MD) by nlin."
    )
    parser.add_argument("--md_csv", type=str, default=DEFAULT_MD_COUNTS_CSV)
    parser.add_argument("--fit_json", type=str, default=DEFAULT_FIT_JSON)
    parser.add_argument("--nlins", type=str, default="64,96,128,160")
    parser.add_argument("--nu", type=float, default=0.5)
    parser.add_argument("--n_trials", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--out_dir", type=str, default=RESULTS_DIR)
    args = parser.parse_args()

    nlins = _parse_nlins(args.nlins)
    md_df = pd.read_csv(args.md_csv)

    required_cols = {"nlin", "linear_length", "count"}
    if not required_cols.issubset(set(md_df.columns)):
        raise ValueError(f"MD CSV must include columns: {required_cols}")

    with open(args.fit_json, "r", encoding="utf-8") as handle:
        fit_data = json.load(handle)
    k1 = float(fit_data["k1"])
    k2 = float(fit_data["k2"])

    # Mapping from MD matrix: nlin 64/96/128/160 -> mlin 6/3/2/1.
    mlin_map = {64: 6, 96: 3, 128: 2, 160: 1}

    summary: dict[str, Any] = {
        "k1": k1,
        "k2": k2,
        "nu": float(args.nu),
        "n_trials": int(args.n_trials),
        "comparisons": {},
    }

    for i, nlin in enumerate(nlins):
        md_sub = md_df[md_df["nlin"] == nlin]
        md_pmf = _normalize_counts(md_sub)
        if not md_pmf:
            print(f"[WARN] No MD counts found for nlin={nlin}; skipping")
            continue

        mlin = int(mlin_map.get(nlin, max(1, int(round(384 / nlin)))))
        rng = np.random.default_rng(args.seed + i)
        dsmc_pmf = simulate_dsmc_length_pmf(
            k1=k1,
            k2=k2,
            nlin=nlin,
            mlin=mlin,
            nu=args.nu,
            n_trials=args.n_trials,
            rng=rng,
        )
        js_div = float(_js_divergence(dsmc_pmf, md_pmf))
        out_plot = _plot_overlay(md_pmf, dsmc_pmf, nlin=nlin, out_dir=args.out_dir)

        summary["comparisons"][str(nlin)] = {
            "nlin": int(nlin),
            "mlin": int(mlin),
            "js_div": js_div,
            "md_pmf": {str(k): float(v) for k, v in sorted(md_pmf.items())},
            "dsmc_pmf": {str(k): float(v) for k, v in sorted(dsmc_pmf.items())},
            "plot": out_plot,
        }
        print(f"nlin={nlin:3d} mlin={mlin} JS={js_div:.6e} plot={out_plot}")

    out_json = save_json(
        summary,
        filename="fit_to_md_pmf_comparison_by_nlin.json",
        out_dir=args.out_dir,
    )
    print(f"saved summary: {out_json}")


if __name__ == "__main__":
    main()
