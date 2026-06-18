# LigMC - Julia Simulation + Python Plotting

This directory keeps the Julia implementation for simulation and sweep execution, while retaining Python plotting scripts for post-processing figures.

## Model summary

At each stage:

1. `mlin` linear polymers of length `nlin` are injected.
2. A Gillespie SSA evolves merge and cyclisation events until no linear chains remain.
3. Every cyclisation creates a new ring that can link to existing rings with Poisson intensity

$$
\mu_t = A\,\frac{n_{\mathrm{ring},t}\,\ell_{\mathrm{cyc}}}{V_{\mathrm{box}}}.
$$

4. The ring graph is updated and gelation is detected from the largest connected component fraction.

The stage box length is rescaled to keep reference monomer density fixed:

$$
\phi_{\mathrm{ref}} = \frac{m_{\mathrm{ring}}n_{\mathrm{ring}} + m_{\mathrm{lin}}n_{\mathrm{lin}}}{L_0^3},
\qquad
L_{\mathrm{stage}} = \left(\frac{M_{\mathrm{total}}}{\phi_{\mathrm{ref}}}\right)^{1/3}.
$$

## Repository layout

```
6_LigMC/
├── simulation_jl/
│   ├── Project.toml
│   ├── Manifest.toml
│   ├── PolymerUtils.jl      # c* scaling, kernels, and constants
│   ├── DSMC.jl              # Gillespie engine (merge + cyclisation)
│   ├── Network.jl           # Ring graph and linking process
│   ├── Main.jl              # Single-trial and multitrial orchestration
│   ├── RunSingle.jl         # CLI for one system configuration
│   ├── SweepGelation.jl     # Grid sweep CLI
│   └── ScalingMrNr.jl       # Scaling analysis for mr*nr vs nr
├── plotting/
│   ├── compare_gel_point_time.py
│   ├── compare_links_per_stage.py
│   ├── plot_gelation_phase_diagram.py
│   └── plot_sim_vs_md_by_nlin.py
├── run_systems.sh           # Batch runner for the MD-matched system list
└── README.md
```

## Setup

```bash
cd simulation_jl
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
cd ..
```

## Usage

Run one system:

```bash
julia --project=simulation_jl simulation_jl/RunSingle.jl \
    --L 80 --mring 78 --nring 512 --mlin 6 --nlin 64 \
    --trials 1000 --n_stages 100 --out_dir ./results
```

Run the standard 16-system matrix:

```bash
bash run_systems.sh
```

Run a parameter sweep:

```bash
julia --project=simulation_jl simulation_jl/SweepGelation.jl \
    --nring_min 256 --nring_max 2048 --nring_step 16 \
    --nlin_min 64 --nlin_max 512 --nlin_step 16 \
    --trials 1000 --n_stages 100 --L 80 --resume \
    --out_dir ./results
```

Plot from sweep results (Python):

```bash
python -m plotting.plot_gelation_phase_diagram \
    --sweep_csv ./results/sweep_summary.csv \
    --out_dir ./results
```

## Scaling of mr*nr with nr

With the concentration construction used in `calculate_polymer_numbers`, ring count satisfies

$$
m_r \propto R_g(n_r)^{-3},
\qquad
R_g(n_r) \propto n_r^{1/2}
\Rightarrow
m_r \propto n_r^{-3/2}.
$$

Therefore,

$$
m_r n_r \propto n_r^{-1/2}.
$$

To compute this directly and fit the exponent numerically:

```bash
julia --project=simulation_jl simulation_jl/ScalingMrNr.jl \
    --nr_min 64 --nr_max 4096 --nr_step 16 --L 80 --ccsr 5.0
```

This writes a CSV with `(nr, mr, mr_nr)` and prints the fitted power law exponent.

## Output files

- `results_all.json`: all trial outputs for a system.
- `summary.csv`: one-row summary with mean/std stages to 50% gel fraction.
- `sweep_summary.csv`: one-row summary per system in a sweep.

## Notes for reproducibility

- Trial seeds are deterministic (`42 + 1337*t`).
- Stage indexing is 1-based in Julia.
- `stages_to_half = nothing` means the trial did not cross 50% by `n_stages`.
