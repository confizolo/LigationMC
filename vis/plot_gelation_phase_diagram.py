"""Create a gelation phase diagram matching the legacy 3_analysis implementation."""

from __future__ import annotations

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants & helpers previously imported from the (now-removed) simulation pkg
# ---------------------------------------------------------------------------
RESULTS_DIR = (
    "/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc"
)
FITTED_A_DEFAULT = 0.20927677484111143


def valence_model(
    l_cyc: float, n_total: float, A: float, box_volume: float = 1.0
) -> float:
    """Mean-field valence: A * n_total * l_cyc / box_volume (> 0 args)."""
    if l_cyc > 0 and n_total > 0 and A > 0:
        return A * n_total * l_cyc / box_volume
    return 0.0


def _fit_power_law(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Fit y = A * x^B in log-log space. Returns (A, B)."""
    valid = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if not np.any(valid):
        return np.nan, np.nan
    x_v, y_v = x[valid], y[valid]
    log_x = np.log10(x_v)
    log_y = np.log10(y_v)
    
    if len(log_x) < 2:
        return np.nan, np.nan
        
    p = np.polyfit(log_x, log_y, 1)
    B = p[0]
    A = 10**p[1]
    return A, B


def plot_gelation_phase_diagram(df: pd.DataFrame, out_dir: str, val_A: float = FITTED_A_DEFAULT) -> str:
    """Create 2x2 phase diagram plot from summary dataframe."""
    os.makedirs(out_dir, exist_ok=True)
    
    if df.empty:
        raise ValueError("Empty dataframe provided to plotter")

    # Filter out nans
    df = df.dropna(subset=["nlin", "nring", "mlin", "mean_stages"])
    df = df[df["mean_stages"] > 0]
    
    if df.empty:
        print("[WARN] No systems with successful gelation found. Plot will be empty.")
        return out_dir
    
    # Calculate valence for the top-left panel using fitted model (assuming B=1.0)
    # L defaults to 80 if not present
    L = df["L"].iloc[0] if "L" in df.columns else 80.0
    box_volume = L**3
    df["valence"] = df.apply(lambda row: valence_model(row["nlin"], row["nring"], val_A, 1.0) / box_volume, axis=1)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    # Common maps
    cmap_nlin = "viridis"
    cmap_nring = "plasma"
    unique_nlins = sorted(df["nlin"].unique())
    
    # ---------------------------------------------------------
    # Panel 1: nlin vs mean_stages
    # ---------------------------------------------------------
    ax = axes[0]
    sc3 = ax.scatter(df["nlin"], df["mean_stages"], c=df["nring"], cmap=cmap_nring, edgecolor="k")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$N_{lin}$")
    ax.set_ylabel("Stages to 50% gelation")
    fig.colorbar(sc3, ax=ax, label=r"$N_{ring}$")
    
    # ---------------------------------------------------------
    # Panel 2: nring vs mean_stages
    # ---------------------------------------------------------
    ax = axes[1]
    sc4 = ax.scatter(df["nring"], df["mean_stages"], c=df["nlin"], cmap=cmap_nlin, edgecolor="k", zorder=2)
    
    for nlin in unique_nlins:
        sub = df[df["nlin"] == nlin]
        if len(sub) > 1:
            x_sub = sub["nring"]
            y_sub = sub["mean_stages"]
            A, B = _fit_power_law(x_sub.to_numpy(), y_sub.to_numpy())
            if not np.isnan(B):
                x_fit = np.linspace(x_sub.min(), x_sub.max(), 100)
                y_fit = A * x_fit**B
                ax.plot(x_fit, y_fit, "k-", alpha=0.3, zorder=1)
                
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$N_{ring}$")
    ax.set_ylabel("Stages to 50% gelation")
    fig.colorbar(sc4, ax=ax, label=r"$N_{lin}$")
    
    # Finalize
    for ax in axes.flat:
        ax.grid(alpha=0.2)
        
    fig.tight_layout()
    out_path = os.path.join(out_dir, "gelation_phase_diagram.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot gelation phase diagram from sweep summary.")
    parser.add_argument("--sweep_csv", type=str, required=True, help="Path to sweep_summary.csv")
    parser.add_argument("--out_dir", type=str, default=RESULTS_DIR)
    parser.add_argument("--val_A", type=float, default=FITTED_A_DEFAULT, help="Valence model A parameter")
    args = parser.parse_args()

    df = pd.read_csv(args.sweep_csv)
    out_path = plot_gelation_phase_diagram(df, out_dir=args.out_dir, val_A=args.val_A)
    print(f"Phase diagram saved to {out_path}")


if __name__ == "__main__":
    main()
