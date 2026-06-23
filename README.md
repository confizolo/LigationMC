# LigMC — Julia Simulation + Python Plotting

Particle DSMC Monte Carlo simulation of topological linking and gelation in
ring-linear polymer blends.  Julia handles simulation and sweeps; Python
handles fitting and post-processing plots.

## Model summary

At each stage:

1. `mlin` linear polymers of length `nlin` are injected.
2. A particle DSMC evolves merge and cyclisation events via acceptance-rejection
   sampling until no linear chains remain.  Each step picks between a merge
   attempt (probability `p_ann`) and a cyclisation attempt (`1 − p_ann`), where
   `p_ann` is set by the ratio of majorant merge and cyclisation rates.  The
   true rate is evaluated and the move is accepted with probability
   `rate / majorant`.
3. Every cyclisation creates a new ring that can link to existing rings with
   Bernoulli probability derived from the Poisson mean

$$
\mu = A\,\frac{n_{\mathrm{ring,target}}\,\ell_{\mathrm{cyc}}}{V_{\mathrm{box}}}.
$$

4. The ring graph is updated and gelation is detected from the largest connected
   component fraction.

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
| `k2` | 7 928.46 | Fitted from MD cyclised-length PMF (nlin=64) via particle DSMC; JS div ≈ 6.5 × 10⁻⁴ |
| `A` | 0.2093 | Fitted from MD valence-by-size data; RMSE = 0.133 |
| `ν` | 0.5 | Ideal chain (Rouse) scaling |

See `_smoke_results/fitted_k1_k2.json` and `_smoke_results/fitted_valence_model.json`.

## Repository layout

```
6_LigMC/
├── src/
│   ├── Project.toml            # Julia dependencies
│   ├── Manifest.toml           # Resolved versions
│   ├── PolymerUtils.jl         # c* scaling, kernels, constants
│   ├── DSMC.jl                 # Particle DSMC engine (merge + cyclisation)
│   ├── Network.jl              # Ring graph and Bernoulli linking
│   ├── Simulation.jl           # Multi-stage trial runner
│   ├── RunSingle.jl            # CLI: single system
│   ├── SweepGelation.jl        # CLI: grid sweep
│   └── ScalingMrNr.jl          # CLI: m_r·n_r scaling analysis
├── vis/
│   ├── fit_dsmc.py             # Fit k2 to MD cyclised-length PMF
│   ├── plot_sim_vs_md_by_nlin.py
│   ├── compare_gel_point_time.py
│   ├── compare_links_per_stage.py
│   └── plot_gelation_phase_diagram.py
├── _smoke_results/             # Diagnostic plots + fitted parameter JSONs
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

Run the standard 16-system matrix:

```bash
julia --project=src src/SweepGelation.jl \
    --nring_min 256 --nring_max 1040 --nring_step 256 \
    --nlin_min 64 --nlin_max 176 --nlin_step 32 \
    --trials 1000 --n_stages 100 --L 80 --resume \
    --out_dir ./results
```

Refit k2 from MD data:

```bash
python vis/fit_dsmc.py
```

Plot DSMC vs MD cyclised-length PMFs:

```bash
python vis/plot_sim_vs_md_by_nlin.py \
    --md_csv /path/to/dist_cyclized_linear_length_by_nlin_all_systems.csv \
    --fit_json _smoke_results/fitted_k1_k2.json \
    --nlins 64,96,128,160
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
