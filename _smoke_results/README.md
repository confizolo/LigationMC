# `_smoke_results/` — Diagnostic Outputs

This directory contains diagnostic plots from early validation runs and symlinks
to the authoritative fitted model parameters stored on cmstore.

## Fitted parameter files (symlinks)

| File | Target | Description |
|---|---|---|
| `fitted_k1_k2.json` | `cmstore02/.../fit_dsmc/fitted_k1_k2.json` | DSMC rate constants: `k1=1.0`, `k2=12933.58`, JS div = 3.4e-4 |
| `fitted_valence_model.json` | `cmstore02/.../fit_poisson/fitted_valence_model.json` | Valence model: `A=0.2093`, `B=1.0`, RMSE = 0.133 |

## Diagnostic plots (legacy)

| File | Description |
|---|---|
| `fit_to_md_pmf_comparison.png` | Cyclised-length PMF: simulation vs MD target (nlin=64) |
| `fit_valence_model_parity.png` | Predicted vs observed valence (parity scatter) |
| `gelation_curves.png` | Gelation fraction vs stage for a sample system |
| `degree_distributions.png` | Degree distribution of the ring graph |
| `compare_scatter.png` | Gelation time parity (sim vs MD) |

## Authoritative results location

All production results (sweep data, per-system trials) live on cmstore:
```
/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc/
├── fit_dsmc/           # k1/k2 fit outputs
├── fit_poisson/        # Valence model fit outputs
├── md_runs/            # 16-system trial results (L=80)
├── sweepL80/           # Full parameter sweep (3627 systems × 1000 trials)
└── sweepL200/          # Extended sweep (L=200)
```
