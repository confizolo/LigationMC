6_LigMC — Focused README
========================

Purpose
-------
This README only covers the `6_LigMC` directory: what each file does, how the data flows, and how to run/debug the key steps.

What you can do with 6_LigMC
----------------------------
- Simulate ring linking events and build a network graph.
- Fit the valence parameter `A` from summary CSVs.
- Compare observed links-per-event with model predictions.

Quick setup (Python 3.10)
-------------------------
```bash
conda create -n newbase python=3.10 -y
conda activate newbase
pip install numpy pandas scipy matplotlib networkx tqdm
```

Core files and what they do
---------------------------
- `main_ligmc.py`
  - Main entry point for running LigMC workflows. It wires together the DSMC engine and the network builder.
  - If you want to run a full simulation, start here and trace how events are created and passed to the network.
  - Box size per stage: `--L` sets the reference box length for stage 0. Each stage recomputes $L$ to keep volume fraction constant based on total monomers.

- `dsmc_engine.py`
  - Defines the event type `CyclisationEvent` used across the directory.
  - Important fields (always check these in debugging):
    - `ring_length`: the newly cyclised ring length (this is `nlin`).
    - `linear_length`: the precursor length (may or may not be used).
    - `links_formed` and `linked_ring_ids`: filled in by `NetworkBuilder.process_cyclisation`.

- `network_builder.py`
  - Implements `NetworkBuilder`, the in-memory graph of rings.
  - Key method: `process_cyclisation(event, A, box_volume=1.0)`
    - Adds the new ring to the graph.
    - For each existing ring, computes:
      - $\mu = A \cdot nring_{target} \cdot nlin / V$
      - $p = 1 - e^{-\mu}$
    - Draws a Bernoulli with probability `p` and adds at most one edge per target.
  - Debug tip: if you see no links, print `event.ring_length`, `nring_target`, `A`, and `box_volume`.

- `polymer_utils.py`
  - Utility functions. `valence_model()` is a legacy helper (A * n_total * l_cyc / V). It is not the active basis used for the current fit.

- `fit_valence_model.py`
  - Fits `A` (and optionally `C`) using summary CSVs that contain `nring`, `mring`, and `nlin`-like columns.
  - Output: JSON with fitted `A` and a PNG overlay plot.
  - If results are off, verify that the basis includes `nring` and `mring`.

- `rate_fitting.py`
  - Fits `k1` and `k2` by matching the **average-length trajectories** from MD (`average_length.txt`) to a Smoluchowski forward model using least squares.
  - Inputs: `--data_root` pointing to MD replica folders, plus `--n_molecules`, `--monomer_length`, and `--volume` (note: this volume is for the MD data; it is not stage-dependent here).
  - Outputs: `fitted_rates.json` with mean and std for `k1`, `k2`, and `kappa`.

- `fit_to_md.py`
  - Fits the **ratio** $\kappa = k_2/k_1$ by matching the cyclized-length PMF from SSA to the MD PMF (default target is `nlin=64`).
  - Default behavior: set `k1=1` (reference scale) and fit only `k2` in log-space to minimize JS-divergence.
  - Outputs: `fitted_k1_k2.json` and `fit_to_md_pmf_comparison.png`.

- `compare_links_per_stage.py`
  - Compares observed links-per-event with model predictions on a grouped CSV.
  - Writes a comparison CSV and optional scatter PNG.

- `tests/smoke_network_test.py`
  - Small harness that repeatedly calls `NetworkBuilder.process_cyclisation` and prints average links and degree histogram.

Data flow (short version)
-------------------------
1) A DSMC runner creates `CyclisationEvent` objects (see `dsmc_engine.py`).
2) Each event is passed into `NetworkBuilder.process_cyclisation` (see `network_builder.py`).
3) The resulting network statistics are summarized into CSVs (outside this folder).
4) `fit_valence_model.py` uses those CSVs to estimate `A`.
5) `compare_links_per_stage.py` uses the fitted `A` to compare predicted vs observed links per event.

How `k1` and `k2` are found
---------------------------
There are two supported routes depending on what MD data you want to match:

1) **Average-length trajectories** (full time series)
  - Use `rate_fitting.py` to least-squares fit `k1` and `k2` to the MD `average_length.txt` traces.
  - This produces `fitted_rates.json` with mean/std values.

2) **Cyclized-length PMF** (distribution shape)
  - Use `fit_to_md.py` to match the cyclized-length PMF from SSA to the MD PMF.
  - The PMF shape depends mainly on the ratio $\kappa = k_2/k_1$.
  - By default it fixes `k1=1` and fits `k2` (log-space) to minimize JS-divergence.

How box size is handled across stages
-------------------------------------
- The reference volume fraction is computed once from stage-0 totals:
  - $\phi_{ref} = M_{0} / L_{0}^3$
  - $M_0 = m_{ring} \cdot n_{ring} + m_{lin} \cdot n_{lin}$
- At the start of each stage, after new linears are added, the box volume is updated using:
  - $M = M_{rings} + m_{lin} \cdot n_{lin}$
  - $V = M / \phi_{ref}$
  - $L = V^{1/3}$
- This updated $L$ is used for all cyclisation events in that stage and is stored in the event timeline.

How to run the main pieces
--------------------------
1) Smoke-test `NetworkBuilder`:

```bash
python tests/smoke_network_test.py
```

2) Fit `A` from a summary CSV:

```bash
python 6_LigMC/fit_valence_model.py --input path/to/summary_all_systems_links_by_size.csv --max_linear_size 200 --outdir path/to/outdir
```

3) Compare observed vs predicted links-per-event:

```bash
python 6_LigMC/compare_links_per_stage.py --csv path/to/summary.csv --model path/to/fitted_valence_model.json --box-volume 512000 --out-csv compare_out.csv --out-png compare_scatter.png
```

Common debugging checks
-----------------------
- Confirm `event.ring_length` is set and non-zero (this is `nlin`).
- Confirm the fitted basis is what you intend: `nring * mring * nlin / V`.
- If nothing links, check scale of `A` and `box_volume` (p can be near 0 or 1).
- Use deterministic RNG in `NetworkBuilder(seed=...)` for reproducible runs.
- Check `box_length_per_stage` in results to verify the box is rescaling as expected.

Minimal debugging snippets
--------------------------
```python
print('nring_target', nring_target, 'nlin', l_cyc, 'box_volume', box_volume)
print('mu', mu, 'p', 1 - math.exp(-mu))
```

Notes
-----
- `nlin` is always the newly cyclised ring length and must be passed in `CyclisationEvent.ring_length`.
- If you need deeper walkthroughs of `main_ligmc.py` or the DSMC runner, tell me the exact path and I will annotate it.
