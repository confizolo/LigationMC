"""Compare MD linking data against Theoretical Formula and LigMC Simulation."""

import argparse
import os
import pickle
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def parse_system_name(dirname: str) -> dict:
    match = re.search(r"L([\d\.]+)_mring(\d+)_nring(\d+)_mlin(\d+)_nlin(\d+)", dirname)
    if not match:
        raise ValueError(f"Could not parse directory name: {dirname}")
    return {
        "L": float(match.group(1)),
        "mring": int(match.group(2)),
        "nring": int(match.group(3)),
        "mlin": int(match.group(4)),
        "nlin": int(match.group(5)),
    }

def get_simulated_events(results_path: str, nring: int, mring: int) -> list[dict]:
    if results_path.endswith(".pkl"):
        with open(results_path, "rb") as f:
            data = pickle.load(f)
    elif results_path.endswith(".json"):
        import json
        with open(results_path, "r") as f:
            data = json.load(f)
    else:
        raise ValueError("Unsupported file format")

    if isinstance(data, dict):
        trials = next(iter(data.values()))
    else:
        trials = data

    events = []
    for trial in trials:
        for event in trial.get("event_timeline", []):
            if event.get("stage") == 1 and event.get("event_type") == "cyclisation":
                events.append({
                    "nring": nring,
                    "mring": mring,
                    "l_cyc": int(event.get("ring_length", event.get("cyclised_length", 0))),
                    "links_formed": int(event.get("links_formed", 0))
                })
    return events

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--md-csv", required=True)
    parser.add_argument("--sim-root", required=True)
    parser.add_argument("--model-json", required=True)
    parser.add_argument("--out-png", required=True)
    args = parser.add_argument_group()
    args = parser.parse_args()

    md_df = pd.read_csv(args.md_csv)

    import json
    with open(args.model_json, "r") as f:
        model = json.load(f)
    A = float(model["A"])

    sim_events = []
    if os.path.isdir(args.sim_root):
        for entry in os.listdir(args.sim_root):
            path = os.path.join(args.sim_root, entry)
            if not os.path.isdir(path):
                continue
            
            try:
                params = parse_system_name(entry)
            except ValueError:
                continue

            results_path_pkl = os.path.join(path, "results_all.pkl")
            results_path_json = os.path.join(path, "results_all.json")
            if os.path.exists(results_path_json):
                results_path = results_path_json
            elif os.path.exists(results_path_pkl):
                results_path = results_path_pkl
            else:
                continue

            sim_events.extend(get_simulated_events(results_path, params["nring"], params["mring"]))

    sim_events_df = pd.DataFrame(sim_events)
    sim_agg = sim_events_df.groupby(["nring", "mring", "l_cyc"]).agg(
        sim_mean=("links_formed", "mean"),
        sim_std=("links_formed", lambda x: np.std(x, ddof=1) if len(x) > 1 else 0.0),
        sim_samples=("links_formed", "count")
    ).reset_index()

    md_df = md_df.rename(columns={"linear_size": "l_cyc"})
    merged = pd.merge(md_df, sim_agg, on=["nring", "l_cyc"], how="inner")

    if merged.empty:
        print("Warning: No matching simulation data found for the MD points!")
        return

    # Expected number of targets is exactly the mring used in the simulation
    # L is fixed at 80.0
    L = 80.0
    
    # Predicted mu per target ring = A * nring * l_cyc / L^3
    merged["mu_target"] = A * merged["nring"] * merged["l_cyc"] / (L ** 3)
    merged["predicted_mean"] = merged["mring"] * (1.0 - np.exp(-merged["mu_target"]))

    fig, ax = plt.subplots(figsize=(6.5, 6))

    merged = merged.sort_values(by=["nring", "l_cyc"])
    
    import matplotlib.cm as cm
    unique_nrings = sorted(merged["nring"].unique())
    colors = cm.viridis(np.linspace(0, 0.9, len(unique_nrings)))
    color_map = {nring: color for nring, color in zip(unique_nrings, colors)}

    for nring, group in merged.groupby("nring"):
        color = color_map[nring]
        
        ax.errorbar(group["l_cyc"], group["avg_links_created"], yerr=group["std_links_created"],
                    fmt="-o", color=color, markersize=7, label=f"MD (nring={nring})")
                    
        ax.errorbar(group["l_cyc"], group["sim_mean"], yerr=group["sim_std"],
                    fmt="--s", color=color, markerfacecolor="white", markersize=7, label=f"LigMC (nring={nring})")
                    
        ax.plot(group["l_cyc"], group["predicted_mean"],
                ":^", color=color, markerfacecolor="white", markersize=7, label=f"Theory (nring={nring})")

    ax.set_xlabel(r"Cyclised Ring Length ($\ell_{cyc}$)", fontsize=12)
    ax.set_ylabel("Links Created per Cyclised Ring", fontsize=12)
    ax.set_title("Linking Events in Stage 1", fontsize=14, weight="bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color="black", lw=2, linestyle="-", marker="o", label="MD (ref)"),
        Line2D([0], [0], color="black", lw=2, linestyle="--", marker="s", markerfacecolor="white", label="LigMC (DSMC)"),
        Line2D([0], [0], color="black", lw=2, linestyle=":", marker="^", markerfacecolor="white", label="Poisson Theory"),
    ]
    for nring in unique_nrings:
        custom_lines.append(Line2D([0], [0], color=color_map[nring], lw=2, label=f"nring={nring}"))
    
    ax.legend(handles=custom_lines, fontsize=10, loc="upper left", bbox_to_anchor=(1.05, 1), frameon=False)
    fig.tight_layout()

    plt.savefig(args.out_png, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.out_png}")
    
    csv_out = args.out_png.replace(".png", ".csv")
    out_df = merged[["nring", "l_cyc", "avg_links_created", "std_links_created", "sim_mean", "sim_std", "predicted_mean"]]
    out_df.to_csv(csv_out, index=False)
    print(f"Wrote {csv_out}")

if __name__ == "__main__":
    main()
