# `parameters/` — Fitted Physics Parameters

This directory securely stores the optimally fitted physics parameters required by the `LigMC` topological simulation and visualisation scripts.

Keeping these lightweight JSON files tracked in the repository ensures that all analysis scripts (like `vis/make_all_comparisons.sh`) can access the exact "best-fit" configuration needed to reproduce the validation plots without needing to re-run the computationally expensive numerical optimizations.

## Files

| File | Description | Generating Script |
|---|---|---|
| `fitted_k1_k2.json` | Contains the intrinsic cyclisation rate constant `k2 = 7928.46` (with `k1 = 1.0`), fitted to exactly match the DSMC PMF to the MD cyclised-length distribution. | `vis/fit_dsmc.py` |
| `fitted_valence_model.json` | Contains the Poisson topological linking prefactor `A = 0.2093`, fitted by minimizing the RMSE between the theoretical linking equation and empirical MD threading counts. | `vis/fit_valence.py` |
