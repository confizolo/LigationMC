# LigMC — Julia Simulation + Python Plotting

Particle DSMC Monte Carlo simulation of topological linking and gelation in
ring-linear polymer blends.  Julia handles simulation and sweeps; Python
handles fitting and post-processing plots.

## Model summary

At each stage:

1. `mlin` linear polymers of length `nlin` are injected in a solution of `mring` ring polymers of length `nring`.
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

See `parameters/fitted_k1_k2.json` and `parameters/fitted_valence_model.json`.

### Fitting the Kinetics Model (k1 and k2)

The kinetic rates for merging ($k_1$) and cyclisation ($k_2$) determine the topological distribution of loop lengths generated during the polymerisation phase.

1. **$k_1$ (Reference Scale):** This parameter is fixed to an arbitrary reference constant ($1.0$). Since we are primarily concerned with the final structural topology (the sequence of merge vs. cyclisation events) rather than the absolute physical kinetics time, fixing $k_1$ simply sets the baseline time scale, and all other physics are driven by the *ratio* of the rates.
2. **$k_2$ (Cyclisation Rate):** This parameter is explicitly fitted using an optimization loop around the Particle DSMC engine. For a reference system (e.g., $n_{lin} = 64, m_{lin} = 6$), the DSMC algorithm simulates the complete polymerisation process to generate a Probability Mass Function (PMF) of the resulting cyclised ring lengths. We then computationally minimize the Jensen-Shannon (JS) divergence between the simulated DSMC PMF and the ground truth PMF extracted from the MD data to find the optimal $k_2$.

### Fitting the Linking Model (Valence Model)

The topological linking parameter $A = 0.2093$ is fitted directly against Molecular Dynamics (MD) tracking data. 

1. **MD Ground Truth:** For each cyclisation event in the MD reference simulations, the actual length of the newly cyclised ring ($\ell_{cyc}$) and the number of topological links formed with existing target rings ($n_{ring}$) are extracted.
2. **Theoretical Expectation:** The expected number of links for each event is modelled as a Poisson process governed by the equation $\mu = A\,\frac{n_{\mathrm{target}}\,\ell_{\mathrm{cyc}}}{V_{\mathrm{box}}}$.
3. **RMSE Minimization:** The parameter $A$ is numerically optimized by minimizing the Root Mean Squared Error (RMSE) between the theoretical expected links and the empirically observed MD link counts across all combinations of $n_{ring}$ and $\ell_{cyc}$.

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
│   ├── fit_valence.py          # Fit A to MD valence data
│   ├── plot_sim_vs_md_by_nlin.py
│   ├── plot_valence_md_comparison.py
│   ├── compare_gel_point_time.py
│   ├── compare_links_per_stage.py
│   ├── plot_gelation_phase_diagram.py
│   └── make_all_comparisons.sh
├── parameters/                 # Diagnostic plots + fitted parameter JSONs
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

Refit parameter A from MD valence data:

```bash
python vis/fit_valence.py
```

Plot DSMC vs MD cyclised-length PMFs:

```bash
python vis/plot_sim_vs_md_by_nlin.py \
    --md_csv /path/to/dist_cyclized_linear_length_by_nlin_all_systems.csv \
    --fit_json parameters/fitted_k1_k2.json \
    --nlins 64,96,128,160
```


## Notes for reproducibility

- Trial seeds are deterministic (`42 + t * 1337`).
- Stage indexing is 1-based in Julia.
- `stages_to_half = nothing` means the trial did not cross 50 % by `n_stages`.
