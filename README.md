# LigMC — Julia Simulation + Python Plotting

DSMC-driven Monte Carlo simulation of topological linking and gelation in
ring-linear polymer blends. Julia handles simulation and sweeps; Python
handles post-processing plots.

## Model summary

At each stage:

1. `mlin` linear polymers of length `nlin` are injected.
2. A Gillespie SSA evolves merge and cyclisation events until no linear chains remain.
3. Every cyclisation creates a new ring that can link to existing rings with
   Bernoulli probability derived from the Poisson mean

$$
\mu = A\,\frac{n_{\mathrm{ring,target}}\,\ell_{\mathrm{cyc}}}{V_{\mathrm{box}}}.
$$

4. The ring graph is updated and gelation is detected from the largest connected component fraction.

The stage box length is rescaled to keep reference monomer density fixed:

$$
\phi_{\mathrm{ref}} = \frac{m_{\mathrm{ring}}n_{\mathrm{ring}} + m_{\mathrm{lin}}n_{\mathrm{lin}}}{L_0^3},
\qquad
L_{\mathrm{stage}} = \left(\frac{M_{\mathrm{total}}}{\phi_{\mathrm{ref}}}\right)^{1/3}.
$$

## Fitted parameters

| Parameter | Value | Provenance |
|---|---|---|
| `k1` | 1.0 | Reference scale (arbitrary) |
| `k2` | 12 933.58 | Fitted from MD cyclised-length PMF (nlin=64); JS div = 3.4 × 10⁻⁴ |
| `A` | 0.2093 | Fitted from MD valence-by-size data; RMSE = 0.133 |
| `ν` | 0.5 | Ideal chain (Rouse) scaling |

See `_smoke_results/fitted_k1_k2.json` and `_smoke_results/fitted_valence_model.json` for full outputs.

## Repository layout

```
6_LigMC/
├── src/
│   ├── Project.toml            # Julia dependencies
│   ├── Manifest.toml           # Resolved versions
│   ├── PolymerUtils.jl         # c* scaling, kernels, constants
│   ├── DSMC.jl                 # Gillespie SSA engine (merge + cyclisation)
│   ├── Network.jl              # Ring graph and Bernoulli linking
│   ├── Simulation.jl           # Multi-stage trial runner
│   ├── RunSingle.jl            # CLI: single system
│   ├── SweepGelation.jl        # CLI: grid sweep
│   └── ScalingMrNr.jl          # CLI: m_r·n_r scaling analysis
├── vis/
│   ├── compare_gel_point_time.py
│   ├── compare_links_per_stage.py
│   ├── plot_gelation_phase_diagram.py
│   └── plot_sim_vs_md_by_nlin.py
├── _smoke_results/             # Diagnostic plots + symlinks to cmstore fits
├── gel_time_compare.csv        # 16-system sim-vs-MD gelation times
├── implementation_plan.md      # v2 design document
├── task.md                     # Execution checklist
├── plan.md                     # Cleanup & completion plan
└── README.md
```

## Setup

```bash
cd src
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
cd ..
```

## Usage

Run one system:

```bash
julia --project=src src/RunSingle.jl \
    --mring 78 --nring 512 --mlin 6 --nlin 64 \
    --trials 1000 --n_stages 100 --out_dir ./results
```

Run the standard 16-system matrix with a sweep restricted to the MD-matched grid:

```bash
julia --project=src src/SweepGelation.jl \
    --nring_min 256 --nring_max 1040 --nring_step 256 \
    --nlin_min 64 --nlin_max 176 --nlin_step 32 \
    --trials 1000 --n_stages 100 --L 80 --resume \
    --out_dir ./results
```

Run a fine-grained parameter sweep:

```bash
julia --project=src src/SweepGelation.jl \
    --nring_min 256 --nring_max 1040 --nring_step 16 \
    --nlin_min 64 --nlin_max 512 --nlin_step 16 \
    --trials 1000 --n_stages 100 --L 80 --resume \
    --out_dir /storage/cmstore02/.../sweepL80
```

Scaling analysis:

```bash
julia --project=src src/ScalingMrNr.jl \
    --nr_min 64 --nr_max 4096 --nr_step 16 --L 80 --ccsr 5.0
```

Plot from sweep results (Python):

```bash
python vis/plot_gelation_phase_diagram.py \
    --sweep_csv /path/to/sweep_summary.csv

python vis/plot_sim_vs_md_by_nlin.py \
    --md_csv /path/to/dist_cyclized_linear_length_by_nlin_all_systems.csv \
    --fit_json _smoke_results/fitted_k1_k2.json
```

## Production results (cmstore)

All production results live on cmstore:

```
/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc/
├── fit_dsmc/           # k1/k2 fit: fitted_k1_k2.json + PMF comparison plots
├── fit_poisson/        # Valence model fit: fitted_valence_model.json + parity plots
├── md_runs/            # 16-system trial results (L=80, 1000 trials each)
├── sweepL80/           # Full parameter sweep (3627 systems × 1000 trials)
└── sweepL200/          # Extended sweep (L=200)
```

## Notes for reproducibility

- Trial seeds are deterministic (`42 + t * 1337`).
- Stage indexing is 1-based in Julia.
- `stages_to_half = nothing` means the trial did not cross 50 % by `n_stages`.
