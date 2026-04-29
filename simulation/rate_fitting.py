"""Fit Smoluchowski rates (k1, k2) from average-length trajectories."""

from __future__ import annotations

import argparse
import glob
import os
from typing import Any

import numpy as np

from .analysis import save_json
from .polymer_utils import RESULTS_DIR

try:
    from scipy.optimize import least_squares
except ImportError as exc:  # pragma: no cover
    raise ImportError("scipy is required for rate fitting: pip install scipy") from exc


def smoluchowski_forward(
    k1: float,
    k2: float,
    n_molecules: int,
    monomer_length: int,
    dt: float,
    n_steps: int,
    alpha: float = 1.0,
    nu: float = 0.6,
    volume: float = 1.0,
) -> np.ndarray:
    """Explicit Euler integration of linear/ring Smoluchowski equations."""
    lengths = monomer_length * (np.arange(n_molecules, dtype=float) + 1.0)
    n_l = np.zeros(n_molecules, dtype=float)
    n_l[0] = n_molecules / volume
    n_r = np.zeros(n_molecules, dtype=float)

    lav_total = np.empty(n_steps, dtype=float)
    lav_total[0] = monomer_length

    kernel_base = np.empty((n_molecules, n_molecules), dtype=float)
    for i in range(n_molecules):
        for j in range(n_molecules):
            ii = i + 1
            jj = j + 1
            kernel_base[i, j] = (ii ** (-alpha) + jj ** (-alpha)) * (ii**nu + jj**nu)

    cyc_scale = np.array([(k + 1) ** (-4.0 * nu) for k in range(n_molecules)], dtype=float)

    for step in range(1, n_steps):
        n_l_new = n_l.copy()
        n_r_new = n_r.copy()

        for k in range(n_molecules):
            for i in range(n_molecules):
                if i + k < n_molecules:
                    n_l_new[k] -= dt * n_l[i] * n_l[k] * k1 * kernel_base[i, k]

                for j in range(n_molecules):
                    if i + j == k:
                        n_l_new[k] += dt * 0.5 * n_l[i] * n_l[j] * k1 * kernel_base[i, j]

            ring_gain = dt * k2 * cyc_scale[k] * n_l[k]
            n_r_new[k] += ring_gain
            n_l_new[k] -= ring_gain

        n_l_new = np.maximum(n_l_new, 0.0)
        n_r_new = np.maximum(n_r_new, 0.0)

        numer = float(np.sum((n_l_new + n_r_new) * lengths))
        denom = float(np.sum(n_l_new + n_r_new))
        lav_total[step] = numer / denom if denom > 0 else lav_total[step - 1]

        n_l = n_l_new
        n_r = n_r_new

    return lav_total


def _load_average_length_files(data_root: str, n_frames: int) -> tuple[np.ndarray, list[np.ndarray]]:
    files = sorted(glob.glob(os.path.join(data_root, "data*", "output", "average_length.txt")))
    if not files:
        files = sorted(glob.glob(os.path.join(data_root, "**", "average_length.txt"), recursive=True))

    traces: list[np.ndarray] = []
    t_ref: np.ndarray | None = None
    for path in files:
        arr = np.loadtxt(path)
        if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] < n_frames:
            continue
        if t_ref is None:
            t_ref = arr[:n_frames, 0]
        traces.append(arr[:n_frames, 1])

    if t_ref is None or not traces:
        raise FileNotFoundError("No valid average_length.txt trajectories were found.")

    return t_ref, traces


def fit_rates(
    data_dirs: list[str],
    n_molecules: int = 200,
    monomer_length: int = 174,
    volume: float = 346.938**3,
    n_frames: int = 1001,
    replicas_per_group: int = 10,
    alpha: float = 1.0,
    nu: float = 0.6,
) -> dict[str, Any]:
    """Fit k1, k2 and kappa from topology-reconstruction average-length outputs."""
    all_rates: list[tuple[float, float]] = []

    for data_root in data_dirs:
        t, traces = _load_average_length_files(data_root, n_frames)
        dt = float(t[1] - t[0])

        for start in range(0, len(traces), replicas_per_group):
            group = traces[start : start + replicas_per_group]
            if not group:
                continue
            mean_trace = np.mean(np.vstack(group), axis=0)

            def residual(k: np.ndarray) -> np.ndarray:
                model = smoluchowski_forward(
                    k1=float(k[0]),
                    k2=float(k[1]),
                    n_molecules=n_molecules,
                    monomer_length=monomer_length,
                    dt=dt,
                    n_steps=n_frames,
                    alpha=alpha,
                    nu=nu,
                    volume=volume,
                )
                return model - mean_trace

            fit = least_squares(
                residual,
                x0=np.array([1e-6, 1e-8], dtype=float),
                bounds=(np.array([0.0, 0.0]), np.array([40.0, 2.0])),
            )
            all_rates.append((float(fit.x[0]), float(fit.x[1])))

    if not all_rates:
        raise RuntimeError("No rate estimates were produced during fitting.")

    rates = np.asarray(all_rates, dtype=float)
    k1_values = rates[:, 0]
    k2_values = rates[:, 1]

    n_density = n_molecules / volume
    kappa_values = 2.0 * k2_values / (n_density * k1_values)

    return {
        "k1": float(np.mean(k1_values)),
        "k2": float(np.mean(k2_values)),
        "kappa": float(np.mean(kappa_values)),
        "k1_std": float(np.std(k1_values, ddof=1)) if len(k1_values) > 1 else 0.0,
        "k2_std": float(np.std(k2_values, ddof=1)) if len(k2_values) > 1 else 0.0,
        "kappa_std": float(np.std(kappa_values, ddof=1)) if len(kappa_values) > 1 else 0.0,
        "n_groups": int(len(all_rates)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit k1/k2 Smoluchowski rates from MD-derived trajectories.")
    parser.add_argument("--data_root", type=str, required=True, help="Root path containing average_length.txt replicas")
    parser.add_argument("--n_molecules", type=int, default=200)
    parser.add_argument("--monomer_length", type=int, default=174)
    parser.add_argument("--volume", type=float, default=346.938**3)
    parser.add_argument("--n_frames", type=int, default=1001)
    parser.add_argument("--replicas_per_group", type=int, default=10)
    parser.add_argument("--out_dir", type=str, default=RESULTS_DIR)
    args = parser.parse_args()

    result = fit_rates(
        data_dirs=[args.data_root],
        n_molecules=args.n_molecules,
        monomer_length=args.monomer_length,
        volume=args.volume,
        n_frames=args.n_frames,
        replicas_per_group=args.replicas_per_group,
    )

    os.makedirs(args.out_dir, exist_ok=True)
    out = save_json(result, filename="fitted_rates.json", out_dir=args.out_dir)

    print("Fitted rates:")
    print(f"k1 = {result['k1']:.6e} +/- {result['k1_std']:.6e}")
    print(f"k2 = {result['k2']:.6e} +/- {result['k2_std']:.6e}")
    print(f"kappa = {result['kappa']:.6e} +/- {result['kappa_std']:.6e}")
    print(f"saved: {out}")


if __name__ == "__main__":
    main()
