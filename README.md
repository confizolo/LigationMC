# 6_LigMC

LigMC contains a staged simulation workflow for ring-link network growth driven by a DSMC gelation model of reactive linear polymers.

## Files

- `main_ligmc.py`: Entry point for staged simulations (CLI + optional multiprocessing).
- `dsmc_engine.py`: Gillespie SSA engine for merge/cyclisation events.
- `network_builder.py`: Graph growth logic for ring-link topology.
- `polymer_utils.py`: Shared physics utilities and output path constants.
- `analysis.py`: Result serialization and plotting helpers.
- `rate_fitting.py`: Standalone least-squares fitting for `k1`, `k2`, and `kappa` from MD-derived `average_length.txt` data.
- `fit_to_md.py`: Fits `k1/k2` by matching cyclized-length PMF to MD target distribution (`nlin=64`).
- `fit_valence_model.py`: Fits concentration-driven linking model `lambda = A * l_cyc^B * N_total` from MD valence data, optimizing only `A` while keeping `B` fixed.
- `simulate_network_growth_progressive_legacy.py`: Legacy model kept for reference.
- `implementation_plan.md`: Refactor plan and design notes.
- `PLAN.md`: High-level implementation context.
- `DSMC/`: Original MATLAB/Python fitting and topology reconstruction utilities.

## Typical Usage

Run a quick staged simulation:

```bash
python main_ligmc.py --mring 50 --nring 256 --mlin 2 --nlin 128 --n_stages 100 --trials 10 --nu 0.5 --val_A 7.393464
```

Note: `main_ligmc.py` defaults now use PMF-fitted rates (`k1=1.0`, `k2=12840.849325710129`) unless you override `--k1/--k2`.

Fit `k1` and `k2` from MD trajectories:

```bash
python rate_fitting.py --data_root /path/to/replicas --n_molecules 200 --monomer_length 174 --volume 41781609.5
```

Fit the PMF-constrained `k1/k2` ratio with Gillespie SSA:

```bash
python fit_to_md.py --nlin 64 --mlin 6 --nu 0.5 --n_trials 5000
```

Fit a single global valence prefactor from MD summary CSV (with B hardcoded to 1):

```bash
python fit_valence_model.py --csv /path/to/summary_all_systems_links_by_size.csv --max_linear_size 256
```

## Outputs

By default, outputs are written to:

`/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc`

Generated artifacts include per-system trial pickles, consolidated results, event timelines, and plots (`gelation_curves.png`, `degree_distributions.png`, `stages_vs_nlin_nring.png`).
