# `src/` — Julia Simulation Modules

This directory contains the Julia implementation of the LigMC model: a particle
DSMC Monte Carlo simulation of topological linking in ring-linear polymer blends.

## Files

| File | Description |
|---|---|
| `PolymerUtils.jl` | Shared constants (`RESULTS_DIR`, reference `Rg` values) and helper functions: `calculate_polymer_numbers`, `smoluchowski_kernel`, `cyclisation_rate`, `valence_model` |
| `DSMC.jl` | Particle DSMC engine.  Defines `Event`, `MergeEvent`, `CyclisationEvent`, `run_dsmc!`.  Uses acceptance-rejection sampling with majorant rates for merge and cyclisation |
| `Network.jl` | Ring graph builder.  Adds rings as nodes, processes cyclisation events with Bernoulli linking, computes largest-component fraction and degree distributions |
| `Simulation.jl` | Multi-stage trial runner.  Manages box-volume rescaling, calls DSMC + Network per stage, tracks gelation metrics |
| `SweepGelation.jl` | CLI entry point for parameter sweeps over `(nring, nlin)` grids.  Saves `sweep_summary.csv` and per-system `results_all.json` |
| `RunSingle.jl` | CLI entry point for a single system configuration.  Saves `results_all.json` and `summary.csv` |
| `ScalingMrNr.jl` | Scaling analysis: computes `m_r · n_r` vs `n_r` and fits a power-law exponent |
| `Project.toml` | Julia package dependencies (ArgParse, Distributions, Graphs, JSON3, ProgressMeter) |
| `Manifest.toml` | Resolved dependency versions |

## Fitted defaults (hard-coded in `Simulation.jl`)

| Parameter | Value | Source |
|---|---|---|
| `k1` | 1.0 | Reference scale |
| `k2` | 7 928.46 | Fitted from cyclised-length PMF via particle DSMC (`parameters/fitted_k1_k2.json`) |
| `A` | 0.2093 | Fitted from MD valence data (`parameters/fitted_valence_model.json`) |
| `ν` | 0.5 | Ideal chain / Rouse scaling |
