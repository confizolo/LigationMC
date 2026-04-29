"""Fit valence models with B hardcoded to 1.

Model used for fitting:
    lambda = A * ((nring * mring * l_cyc) / V_box)

Optional offset form:
    lambda = A * ((nring * mring * l_cyc) / V_box) + C
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any

# Ensure repo root is on sys.path for direct script execution.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import matplotlib.pyplot as plt
import numpy as np

from simulation.analysis import save_json
from simulation.polymer_utils import RESULTS_DIR

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover
    raise ImportError("pandas is required for fit_valence_model.py: pip install pandas") from exc

DEFAULT_MD_CSV = (
    "/storage/cmstore02/groups/TAPLab/fconforto-projects/"
    "fconforto-olympic-gels/results/histories/summary_all_systems_links_by_size.csv"
)


def _resolve_column(df: pd.DataFrame, candidates: list[str]) -> str:
    lower_map = {col.lower(): col for col in df.columns}
    for cand in candidates:
        key = cand.lower()
        if key in lower_map:
            return lower_map[key]
    raise KeyError(f"Could not find any of columns: {candidates}")


def _infer_n_total_from_companion(md_data: pd.DataFrame, csv_path: str | None = None) -> pd.Series | None:
    """Infer ring-count n_total from companion summary metadata when mring is absent.

    The grouped links CSV commonly keeps `nring` (ring polymerization index) but drops
    `mring` (number of rings). The valence model needs the latter for density.
    """
    if "nring" not in md_data.columns:
        return None

    if not csv_path:
        return None

    companion = os.path.join(os.path.dirname(os.path.abspath(csv_path)), "summary_rep.csv")
    if not os.path.exists(companion):
        return None

    rep = pd.read_csv(companion)
    if "nring" not in rep.columns or "mring" not in rep.columns:
        return None

    pairs = rep[["nring", "mring"]].dropna().drop_duplicates()
    if pairs.empty:
        return None

    # Require one unique mring per nring to avoid silently mixing inconsistent metadata.
    counts = pairs.groupby("nring")["mring"].nunique()
    ambiguous = counts[counts > 1]
    if not ambiguous.empty:
        bad = ", ".join(str(int(v)) for v in sorted(ambiguous.index.astype(int).tolist()))
        raise ValueError(
            "Cannot infer n_total from companion summary_rep.csv because nring->mring "
            f"mapping is ambiguous for nring values: {bad}"
        )

    mapping = pairs.groupby("nring")["mring"].first().to_dict()
    inferred = pd.to_numeric(md_data["nring"], errors="coerce").map(mapping)
    return pd.to_numeric(inferred, errors="coerce")


def _resolve_columns(md_data: pd.DataFrame) -> dict[str, str | None]:
    """Map semantic column names to actual column names in the DataFrame."""
    linear_col = _resolve_column(md_data, ["ring_size", "linear_size", "nlin", "length", "l_cyc"])
    links_col = _resolve_column(md_data, ["avg_links_created", "avg_links", "mean_links", "lambda"])

    n_total_col = None
    try:
        n_total_col = _resolve_column(md_data, ["mring", "n_total", "ring_count"])
    except KeyError:
        pass

    n_samples_col = None
    try:
        n_samples_col = _resolve_column(md_data, ["n_samples", "samples", "count", "n"])
    except KeyError:
        pass

    std_col = None
    try:
        std_col = _resolve_column(md_data, ["std_links_created", "std_links", "std", "sigma"])
    except KeyError:
        pass

    return {
        "linear_size": linear_col,
        "avg_links_created": links_col,
        "n_total": n_total_col,
        "n_samples": n_samples_col,
        "std_links_created": std_col,
    }


def _infer_missing_columns(df: pd.DataFrame, md_data: pd.DataFrame, csv_path: str | None) -> pd.DataFrame:
    """Fill n_total from companion metadata if absent."""
    if "n_total" not in df.columns:
        inferred_n_total = _infer_n_total_from_companion(md_data, csv_path=csv_path)
        if inferred_n_total is not None:
            df["n_total"] = inferred_n_total
        elif "nring" in df.columns:
            # Last resort: keep backward-compatible behavior but make the fallback explicit.
            df["n_total"] = pd.to_numeric(df["nring"], errors="coerce")
            print(
                "[WARN] Falling back to nring as n_total because mring/n_total is unavailable. "
                "This can invert expected slope ordering if nring is ring size rather than ring count."
            )

    # Keep nring for overlay split and basis construction.
    lower_cols = {c.lower(): c for c in md_data.columns}
    if "nring" in lower_cols:
        df["nring"] = pd.to_numeric(md_data[lower_cols["nring"]], errors="coerce")
    elif "nring" in df.columns:
        df["nring"] = pd.to_numeric(df["nring"], errors="coerce")

    if "nring" not in df.columns:
        raise KeyError(
            "Could not find nring column required by basis: (nring * mring * l_cyc) / V_box"
        )
    return df


def _filter_and_coerce(df: pd.DataFrame, max_linear_size: int) -> pd.DataFrame:
    """Coerce types, drop invalid rows, and apply size filter."""
    df["linear_size"] = pd.to_numeric(df["linear_size"], errors="coerce")
    df["avg_links_created"] = pd.to_numeric(df["avg_links_created"], errors="coerce")
    df["n_total"] = pd.to_numeric(df["n_total"], errors="coerce")
    
    if "n_samples" not in df.columns:
        df["n_samples"] = 1.0
    else:
        df["n_samples"] = pd.to_numeric(df["n_samples"], errors="coerce").fillna(1.0)
        
    if "std_links_created" in df.columns:
        df["std_links_created"] = pd.to_numeric(df["std_links_created"], errors="coerce")

    df = df.dropna(subset=["linear_size", "avg_links_created", "n_total", "nring"])
    df = df[df["linear_size"] > 0]
    df = df[df["n_total"] > 0]
    df = df[df["nring"] > 0]
    df = df[df["linear_size"] <= float(max_linear_size)]

    if df.empty:
        raise ValueError("No valid rows remain after filtering valence data.")

    return df.reset_index(drop=True)


def _prepare_dataframe(md_data: pd.DataFrame, max_linear_size: int = 256, csv_path: str | None = None) -> pd.DataFrame:
    col_map = _resolve_columns(md_data)
    
    keep_cols = [col_map["linear_size"], col_map["avg_links_created"]]
    if col_map["n_total"]:
        keep_cols.append(col_map["n_total"])
    elif "n_total" in md_data.columns:
        keep_cols.append("n_total")
    elif "nring" in md_data.columns:
        keep_cols.append("nring")
    else:
        raise KeyError(
            "Could not find a density/count column. Expected one of: mring, n_total, ring_count "
            "(or nring with companion summary_rep.csv metadata)."
        )
        
    if col_map["n_samples"]:
        keep_cols.append(col_map["n_samples"])
    if col_map["std_links_created"]:
        keep_cols.append(col_map["std_links_created"])

    rename_map = {
        col_map["linear_size"]: "linear_size",
        col_map["avg_links_created"]: "avg_links_created",
    }
    if col_map["n_total"]: rename_map[col_map["n_total"]] = "n_total"
    if col_map["n_samples"]: rename_map[col_map["n_samples"]] = "n_samples"
    if col_map["std_links_created"]: rename_map[col_map["std_links_created"]] = "std_links_created"

    df = md_data[keep_cols].rename(columns=rename_map)
    df = _infer_missing_columns(df, md_data, csv_path)
    return _filter_and_coerce(df, max_linear_size)


def _basis_from_df(df: pd.DataFrame, box_volume: float) -> np.ndarray:
    if box_volume <= 0.0:
        raise ValueError("box_volume must be > 0")
    linear_size = df["linear_size"].to_numpy(dtype=float)
    n_total = df["n_total"].to_numpy(dtype=float)
    nring = df["nring"].to_numpy(dtype=float)
    return (nring * n_total * linear_size) / float(box_volume)


def fit_global_A(
    md_data: pd.DataFrame,
    box_volume: float = 80.0**3,
    max_linear_size: int = 256,
    csv_path: str | None = None,
) -> dict[str, Any]:
    df = _prepare_dataframe(md_data, max_linear_size=max_linear_size, csv_path=csv_path)

    observed = df["avg_links_created"].to_numpy(dtype=float)
    n_samples = np.maximum(df["n_samples"].to_numpy(dtype=float), 1.0)
    std_links = (
        pd.to_numeric(df["std_links_created"], errors="coerce").to_numpy(dtype=float)
        if "std_links_created" in df.columns
        else np.full_like(observed, np.nan, dtype=float)
    )

    # Prefer empirical SEM when present; fallback to Poisson-inspired uncertainty.
    sem = std_links / np.sqrt(n_samples)
    poisson_sem = np.sqrt(np.maximum(observed, 1e-9) / n_samples)
    sigma = np.where(np.isfinite(sem) & (sem > 0.0), sem, poisson_sem)
    sigma = np.maximum(sigma, 1e-9)

    basis = _basis_from_df(df, box_volume=box_volume)
    weights = 1.0 / (sigma * sigma)

    # Weighted least-squares closed form for the only free coefficient A.
    denom = float(np.sum(weights * basis * basis))
    if denom <= 0.0:
        raise ValueError("Degenerate weighted least-squares system while fitting A.")
    A_fit = float(np.sum(weights * basis * observed) / denom)
    A_fit = max(0.0, A_fit)

    pred_fit = A_fit * basis
    rmse = float(np.sqrt(np.mean((pred_fit - observed) ** 2)))

    return {"A": A_fit, "rmse": rmse, "df": df}


def fit_global_A_plus_C(
    md_data: pd.DataFrame,
    box_volume: float = 80.0**3,
    max_linear_size: int = 256,
    csv_path: str | None = None,
) -> dict[str, Any]:
    """Fit lambda = A * basis + C with weighted least squares."""
    df = _prepare_dataframe(md_data, max_linear_size=max_linear_size, csv_path=csv_path)

    observed = df["avg_links_created"].to_numpy(dtype=float)
    n_samples = np.maximum(df["n_samples"].to_numpy(dtype=float), 1.0)
    std_links = (
        pd.to_numeric(df["std_links_created"], errors="coerce").to_numpy(dtype=float)
        if "std_links_created" in df.columns
        else np.full_like(observed, np.nan, dtype=float)
    )

    sem = std_links / np.sqrt(n_samples)
    poisson_sem = np.sqrt(np.maximum(observed, 1e-9) / n_samples)
    sigma = np.where(np.isfinite(sem) & (sem > 0.0), sem, poisson_sem)
    sigma = np.maximum(sigma, 1e-9)

    basis = _basis_from_df(df, box_volume=box_volume)
    w = 1.0 / (sigma * sigma)

    # Weighted normal equations for [A, C].
    s_bb = float(np.sum(w * basis * basis))
    s_b1 = float(np.sum(w * basis))
    s_11 = float(np.sum(w))
    t_b = float(np.sum(w * basis * observed))
    t_1 = float(np.sum(w * observed))

    mat = np.array([[s_bb, s_b1], [s_b1, s_11]], dtype=float)
    rhs = np.array([t_b, t_1], dtype=float)
    det = float(np.linalg.det(mat))
    if abs(det) <= 1e-18:
        raise ValueError("Degenerate weighted least-squares system while fitting A and C.")
    A_fit, C_fit = np.linalg.solve(mat, rhs)

    pred_fit = (A_fit * basis) + C_fit
    rmse = float(np.sqrt(np.mean((pred_fit - observed) ** 2)))

    return {"A": float(A_fit), "C": float(C_fit), "rmse": rmse, "df": df}


def _plot_overlay_by_nring(df: pd.DataFrame, A: float, box_volume: float, out_dir: str, C: float = 0.0, suffix: str = "") -> str:
    if "nring" not in df.columns:
        raise ValueError("Cannot make nring overlay plot because 'nring' column is not available.")

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    cmap = plt.get_cmap("tab10")

    for i, nring in enumerate(sorted(df["nring"].dropna().unique().astype(int))):
        sub = df[df["nring"] == nring].sort_values("linear_size")
        if sub.empty:
            continue

        x = sub["linear_size"].to_numpy(dtype=float)
        y = sub["avg_links_created"].to_numpy(dtype=float)
        n_samples = np.maximum(sub["n_samples"].to_numpy(dtype=float), 1.0)
        std_links = (
            pd.to_numeric(sub["std_links_created"], errors="coerce").to_numpy(dtype=float)
            if "std_links_created" in sub.columns
            else np.full_like(y, np.nan, dtype=float)
        )
        sem = std_links / np.sqrt(n_samples)
        poisson_sem = np.sqrt(np.maximum(y, 1e-9) / n_samples)
        yerr = np.where(np.isfinite(sem) & (sem > 0.0), sem, poisson_sem)

        color = cmap(i % 10)
        ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=3, markersize=5, color=color, label=f"MD nring={nring}")

        n_total_curve = float(np.nanmedian(sub["n_total"].to_numpy(dtype=float)))
        n_total_curve = max(n_total_curve, 1.0)

        xfit = np.linspace(float(np.min(x)), float(np.max(x)), 200)
        basis_fit = (float(nring) * n_total_curve * xfit) / float(box_volume)
        yfit = C + (A * basis_fit)
        ax.plot(
            xfit,
            yfit,
            "--",
            color=color,
            alpha=0.95,
            label=f"Model nring={nring} (mring={n_total_curve:.1f}, A={A:.3e}, C={C:.3e})",
        )

    ax.set_xlabel("Cyclised Ring Size (monomers)")
    ax.set_ylabel("Average Links Created")
    ax.set_title("MD Points vs Valence Model (B=1)")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    filename = f"fit_valence_overlay_by_nring{suffix}.png" if suffix else "fit_valence_overlay_by_nring.png"
    out_path = os.path.join(out_dir, filename)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit global A (and optional C) for basis (nring*mring*l_cyc)/V.")
    parser.add_argument("--csv", type=str, default=DEFAULT_MD_CSV, help="Path to MD summary CSV")
    parser.add_argument("--max_linear_size", type=int, default=256)
    parser.add_argument("--box_length", type=float, default=80.0, help="Box side length used for V_box = box_length^3")
    parser.add_argument("--out_dir", type=str, default=RESULTS_DIR)
    parser.add_argument(
        "--fit_offset",
        action="store_true",
        help="Also fit lambda = A * ((nring*mring*l_cyc)/V_box) + C and report improvement.",
    )
    args = parser.parse_args()

    df_raw = pd.read_csv(args.csv)
    box_volume = float(args.box_length**3)
    fit = fit_global_A(
        df_raw,
        box_volume=box_volume,
        max_linear_size=int(args.max_linear_size),
        csv_path=args.csv,
    )
    df = fit["df"]

    fit_offset: dict[str, Any] | None = None
    if args.fit_offset:
        fit_offset = fit_global_A_plus_C(
            df_raw,
            box_volume=box_volume,
            max_linear_size=int(args.max_linear_size),
            csv_path=args.csv,
        )

    payload = {
        "A": float(fit["A"]),
        "B": 1.0,
        "rmse": float(fit["rmse"]),
        "n_points": int(len(df)),
        "csv": args.csv,
        "box_length": float(args.box_length),
        "box_volume": box_volume,
        "max_linear_size": int(args.max_linear_size),
    }
    if fit_offset is not None:
        payload["A_with_offset"] = float(fit_offset["A"])
        payload["C_with_offset"] = float(fit_offset["C"])
        payload["rmse_with_offset"] = float(fit_offset["rmse"])

    os.makedirs(args.out_dir, exist_ok=True)
    out_json = save_json(payload, filename="fitted_valence_model.json", out_dir=args.out_dir)
    out_overlay = _plot_overlay_by_nring(df, A=float(fit["A"]), box_volume=box_volume, out_dir=args.out_dir)
    out_overlay_offset = None
    if fit_offset is not None:
        out_overlay_offset = _plot_overlay_by_nring(
            fit_offset["df"],
            A=float(fit_offset["A"]),
            C=float(fit_offset["C"]),
            box_volume=box_volume,
            out_dir=args.out_dir,
            suffix="_with_offset",
        )

    print(f"A      = {fit['A']:.6e}")
    print("B      = 1.000000")
    print(f"RMSE   = {fit['rmse']:.6f}")
    if fit_offset is not None:
        print(f"A(+C)  = {fit_offset['A']:.6e}")
        print(f"C      = {fit_offset['C']:.6e}")
        print(f"RMSE+C = {fit_offset['rmse']:.6f}")
    print(f"saved  = {out_json}")
    print(f"plot   = {out_overlay}")
    if out_overlay_offset is not None:
        print(f"plot+C = {out_overlay_offset}")


if __name__ == "__main__":
    main()
