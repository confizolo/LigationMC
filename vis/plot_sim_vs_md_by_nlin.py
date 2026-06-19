"""Create simulation-vs-MD cyclized-length PMF overlays for multiple nlin values."""

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
# Constants & helpers previously imported from the (now-removed) simulation pkg
# ---------------------------------------------------------------------------
RESULTS_DIR = (
    "/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc"
)
FITTED_K1_DEFAULT = 1.0
FITTED_K2_DEFAULT = 12933.579888871815


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


def simulate_length_pmf(
    k1: float,
    k2: float,
    nlin: int,
    mlin: int,
    nu: float,
    n_trials: int,
    rng: np.random.Generator,
) -> dict[int, float]:
    """Gillespie SSA producing the cyclised-ring-length PMF.

    Starts with *mlin* linears of length *nlin*. Merge and cyclisation
    propensities follow the Smoluchowski kernel with exponent *nu* and
    cyclisation exponent -4*nu.

    Returns a normalised PMF dict {length: probability}.
    """
    alpha = 1.0  # fixed exponent in the merge kernel
    counts: dict[int, int] = {}

    for _ in range(n_trials):
        # populations: length -> count of linear chains with that length
        pops: dict[int, int] = {nlin: mlin}

        while True:
            # Collect species list
            species = [(length, cnt) for length, cnt in pops.items() if cnt > 0]
            if not species:
                break

            # --- build propensity list ---
            propensities: list[tuple[float, str, Any]] = []

            # Merge propensities: k1 * ni * nj * K(i,j)
            for idx_a, (la, na) in enumerate(species):
                for idx_b in range(idx_a, len(species)):
                    lb, nb = species[idx_b]
                    if idx_a == idx_b:
                        if na < 2:
                            continue
                        factor = na * (na - 1)  # ordered pairs, matching Julia DSMC
                    else:
                        factor = na * nb
                    if factor <= 0:
                        continue
                    kernel = (la ** (-alpha) + lb ** (-alpha)) * (la**nu + lb**nu)
                    rate = k1 * factor * kernel
                    if rate > 0:
                        propensities.append((rate, "merge", (la, lb)))

            # Cyclisation propensities: ni * k2 * length^(-4*nu)
            for la, na in species:
                if na <= 0:
                    continue
                rate = na * k2 * la ** (-4.0 * nu)
                if rate > 0:
                    propensities.append((rate, "cycle", (la,)))

            if not propensities:
                break

            total_rate = sum(r for r, _, _ in propensities)
            if total_rate <= 0:
                break

            # Gillespie draw
            r1 = rng.random()
            threshold = r1 * total_rate
            cumsum = 0.0
            chosen = propensities[-1]
            for prop in propensities:
                cumsum += prop[0]
                if cumsum >= threshold:
                    chosen = prop
                    break

            _, event_type, payload = chosen

            if event_type == "merge":
                la, lb = payload
                # Remove one of each reactant
                pops[la] -= 1
                if pops[la] <= 0:
                    del pops[la]
                pops[lb] = pops.get(lb, 0) - 1
                if pops[lb] <= 0:
                    del pops[lb]
                # Add merged product
                new_len = la + lb
                pops[new_len] = pops.get(new_len, 0) + 1
            else:  # cycle
                la = payload[0]
                pops[la] -= 1
                if pops[la] <= 0:
                    del pops[la]
                counts[la] = counts.get(la, 0) + 1

    # Normalise to PMF
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


def _plot_overlay(md_pmf: dict[int, float], sim_pmf: dict[int, float], nlin: int, out_dir: str) -> str:
    keys = sorted(set(md_pmf) | set(sim_pmf))
    x = np.arange(len(keys), dtype=float)
    width = 0.42

    y_md = np.array([md_pmf.get(k, 0.0) for k in keys], dtype=float)
    y_sim = np.array([sim_pmf.get(k, 0.0) for k in keys], dtype=float)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(x - width / 2.0, y_md, width=width, label="MD target", alpha=0.85)
    ax.bar(x + width / 2.0, y_sim, width=width, label="LigMC SSA", alpha=0.85)
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
    parser = argparse.ArgumentParser(description="Generate PMF overlays (simulation vs MD) by nlin.")
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

    # Mapping from your MD matrix: nlin 64/96/128/160 -> mlin 6/3/2/1.
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
        sim_pmf = simulate_length_pmf(
            k1=k1,
            k2=k2,
            nlin=nlin,
            mlin=mlin,
            nu=args.nu,
            n_trials=args.n_trials,
            rng=rng,
        )
        js_div = float(_js_divergence(sim_pmf, md_pmf))
        out_plot = _plot_overlay(md_pmf, sim_pmf, nlin=nlin, out_dir=args.out_dir)

        summary["comparisons"][str(nlin)] = {
            "nlin": int(nlin),
            "mlin": int(mlin),
            "js_div": js_div,
            "md_pmf": {str(k): float(v) for k, v in sorted(md_pmf.items())},
            "sim_pmf": {str(k): float(v) for k, v in sorted(sim_pmf.items())},
            "plot": out_plot,
        }
        print(f"nlin={nlin:3d} mlin={mlin} JS={js_div:.6e} plot={out_plot}")

    out_json = save_json(summary, filename="fit_to_md_pmf_comparison_by_nlin.json", out_dir=args.out_dir)
    print(f"saved summary: {out_json}")


if __name__ == "__main__":
    main()
