# LigMC Fitting: Explicit Description

This note documents exactly how simulation parameters are fitted to MD data in the current `6_LigMC` workflow.

## 1) Cyclization PMF fit (`k1`, `k2`)

### Goal
Fit the DSMC rates so the simulated cyclized-length PMF matches the MD PMF.

### MD target used for fitting
The primary fitting target is the `nlin=64` PMF (hard-coded in `fit_to_md.py`):

- 64: 0.7799
- 128: 0.1557
- 192: 0.0371
- 256: 0.0222
- 320: 0.0052
- 384: 0.0015

### Forward model
For candidate `(k1, k2)`:

1. Build one stage with `mlin` linears of length `nlin`.
2. Run Gillespie SSA (`DSMCEngine`) until linears are exhausted.
3. Record all cyclization event ring lengths.
4. Repeat over `n_trials` and form the simulated PMF.

### Objective function
Distance between simulated and MD PMFs is Jensen-Shannon divergence:

`JS(P_sim, P_md)`

The fitter minimizes this divergence.

### Optimization setup
Current default fitting mode is 1D:

- `k1` fixed (`k1_fixed`, default 1.0)
- optimize `k2`
- variable is `log10(k2)`

Implementation details:

- deterministic objective evaluation by re-seeding RNG per candidate point
- objective cache for repeated points
- coarse log-grid scan for robust initialization
- local Nelder-Mead refinement

So effectively:

`k2* = argmin_k2 JS(P_sim(k1_fixed, k2), P_md)`

and ratio is:

`kappa = k2 / k1`

### Outputs
- `fitted_k1_k2.json`
- `fit_to_md_pmf_comparison.png`

## 2) Linking Poisson fit (A-only, B fixed)

### Linking model used in simulation
Current linking intensity per cyclization event is:

`lambda = A * ((N_total * l_cyc) / V_box)^B`

where:

- `N_total`: number of available ring nodes in the pool
- `l_cyc`: cyclized linear length for that event
- `V_box = L^3`
- `A`: fitted prefactor
- `B`: fixed exponent (default 1.0)

### MD data used
`summary_all_systems_links_by_size.csv`, with columns like:

- `nring` (used as `N_total` proxy in fitting)
- `linear_size` (used as `l_cyc`)
- `avg_links_created` (target lambda)
- `std_links_created`, `n_samples` (for uncertainty weighting)

### Objective / estimator
Only `A` is fitted, with `B` fixed.

Let:

`x_i = ((N_i * l_i) / V_box)^B`
`y_i = observed_avg_links_i`

Weights use standard error when available:

- `sigma_i = std_links_created_i / sqrt(n_samples_i)`
- fallback: `sqrt(y_i / n_samples_i)`
- `w_i = 1 / sigma_i^2`

Then `A` is computed by weighted least squares (closed form):

`A* = (sum_i w_i x_i y_i) / (sum_i w_i x_i^2)`

### Reported metric
- RMSE between predicted and observed average links.

### Outputs
- `fitted_valence_model.json`
- `fit_valence_model_parity.png`

## 3) Runtime use in `main_ligmc.py`

- The simulation uses fitted `k1`, `k2` defaults (unless overridden by CLI).
- During each cyclization event, `network_builder.process_cyclisation(...)` computes `lambda` from the model above and draws links from `Poisson(lambda)`.
- `box_volume` is passed as `L^3` from runtime config.

## 4) Evaluation vs fitting (important)

`plot_sim_vs_md_by_nlin.py` does **not** fit parameters. It only:

1. reads fitted `k1/k2`
2. simulates PMFs for selected `nlin`
3. compares to MD PMFs and reports JS

So it is a validation/evaluation script, not an optimizer.

## 5) Practical interpretation of recent behavior

- PMF fit quality improved strongly after stabilizing the objective.
- Linking A-only fit quality remained poor (high RMSE), suggesting model-form mismatch rather than optimizer failure.
- Rewriting concentration as `(N*l)/V` mostly rescales `A` when `V` is fixed; it does not by itself guarantee lower RMSE.
