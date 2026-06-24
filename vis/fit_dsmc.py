"""Fit the DSMC cyclisation rate constant k2 to MD cyclised-length data.

Minimises the Jensen-Shannon divergence between the DSMC-simulated
cyclised-ring-length PMF and the MD target PMF (nlin=64, mlin=6).

Outputs a JSON file with the fitted {k1, k2, js_div}.
"""

from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


def _kl_divergence(p: dict, q: dict) -> float:
    kl = 0.0
    for k in p:
        pk = p[k]
        qk = q.get(k, 0.0)
        if pk > 0.0 and qk > 0.0:
            kl += pk * math.log(pk / qk)
    return kl


def _js_divergence(p: dict, q: dict) -> float:
    all_keys = set(p) | set(q)
    m = {}
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
    """Particle DSMC producing the cyclised-ring-length PMF."""
    alpha = 1.0
    density = k1 * mlin
    counts: dict[int, int] = {}

    for _ in range(n_trials):
        masses = [nlin] * mlin
        n_chains = mlin

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


def get_md_target_pmf(nlin: int = 64) -> dict[int, float]:
    md_csv = (
        "/storage/cmstore02/groups/TAPLab/fconforto-projects/"
        "fconforto-olympic-gels/results/histories/"
        "dist_cyclized_linear_length_by_nlin_all_systems.csv"
    )
    md_df = pd.read_csv(md_csv)
    sub = md_df[md_df["nlin"] == nlin]
    counts = sub.groupby("linear_length", as_index=True)["count"].sum().to_dict()
    total = float(sum(max(v, 0.0) for v in counts.values()))
    return {int(k): float(v / total) for k, v in counts.items() if v > 0}


def main() -> None:
    md_pmf = get_md_target_pmf(nlin=64)
    print("MD Target PMF for nlin=64:", md_pmf)

    k1 = 1.0
    nu = 0.5
    nlin = 64
    mlin = 6
    n_trials = 20000

    def objective(k2: float) -> float:
        if k2 <= 0:
            return 1e6
        rng = np.random.default_rng(42)
        sim_pmf = simulate_dsmc_length_pmf(k1, k2, nlin, mlin, nu, n_trials, rng)
        js = _js_divergence(sim_pmf, md_pmf)
        print(f"k2 = {k2:.2f} -> JS = {js:.6e}")
        return js

    print("\nStarting optimization...")
    res = minimize_scalar(
        objective, bounds=(5000, 15000), method="bounded",
        options={"xatol": 1.0, "maxiter": 50},
    )

    print("\nOptimization Complete:")
    print("Best k2:", res.x)
    print("Best JS:", res.fun)

    out_dict = {"k1": k1, "k2": res.x, "js_div": res.fun, "method": "ParticleDSMC"}

    os.makedirs("./results", exist_ok=True)
    out_path = "./parameters/fitted_k1_k2.json"
    with open(out_path, "w") as f:
        json.dump(out_dict, f, indent=2)
    print(f"Saved fitted parameters to {out_path}")


if __name__ == "__main__":
    main()
