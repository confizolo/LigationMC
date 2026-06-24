# `vis/` — Python Plotting & Fitting Scripts

Post-processing visualisation and fitting scripts for LigMC simulation results.
All scripts are self-contained (no dependency on external simulation packages).

## Files

| File | Description |
|---|---|
| `fit_dsmc.py` | Fits the cyclisation rate constant `k2` by minimising the JS divergence between the DSMC PMF and the MD target (nlin=64).  Outputs `fitted_k1_k2.json` |
| `fit_valence.py` | Fits the Poisson topological linking parameter `A` by minimizing the RMSE against MD link data. Outputs `fitted_valence_model.json` |
| `plot_sim_vs_md_by_nlin.py` | PMF overlay: DSMC vs MD cyclised-length distributions across multiple `nlin` values |
| `compare_gel_point_time.py` | Scatter plot of LigMC gelation time vs MD reference (parity plot with error bars) |
| `compare_links_per_stage.py` | Compares observed links-per-stage with valence-model predictions |
| `plot_valence_md_comparison.py` | Plots stage-1 links created per cyclised ring against the fitted Poisson model |
| `plot_gelation_phase_diagram.py` | 2-panel gelation phase diagram (`nlin` vs stages, `nring` vs stages) from `sweep_summary.csv` |

## Dependencies

- Python 3.10+
- `numpy`, `matplotlib`, `pandas` (standard scientific stack)
- `scipy` (for `fit_dsmc.py` and `plot_sim_vs_md_by_nlin.py`)

## Usage examples

```bash
# Fit k2 to MD cyclised-length data
python vis/fit_dsmc.py

# PMF overlay by nlin
python vis/plot_sim_vs_md_by_nlin.py \
    --md_csv /path/to/dist_cyclized_linear_length_by_nlin_all_systems.csv \
    --fit_json parameters/fitted_k1_k2.json \
    --nlins 64,96,128,160

# Gelation parity plot
python vis/compare_gel_point_time.py \
    --gel-file /path/to/gel_point_time.txt \
    --results-root /path/to/per-system-results \
    --out-csv gel_time_compare.csv \
    --out-png gel_time_compare.png

# Gelation phase diagram from sweep
python vis/plot_gelation_phase_diagram.py \
    --sweep_csv /path/to/sweep_summary.csv
```
