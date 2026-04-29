# LigMC — Ligation Monte Carlo

Staged simulation of ring-polymer network growth driven by reactive linear polymers, using a Gillespie SSA (DSMC) engine with concentration-dependent Poisson linking calibrated to MD data.

## Physical model

A dense melt of ring polymers is threaded by sparse reactive linear polymers. Linears merge end-to-end (Smoluchowski coagulation) and eventually cyclise to form new rings. Each cyclisation event may topologically link the new ring to existing ones, building a network. The simulation tracks this network and identifies the gelation transition (largest connected component ≥ 50%).

### Staged growth

At each stage, `mlin` fresh linears of length `nlin` are injected. The DSMC engine evolves the linear population until all are consumed (cyclised). After each cyclisation, `NetworkBuilder` draws topological links and updates the ring graph. The box length `L` is rescaled per stage to keep the volume fraction constant:

$$\phi_{\text{ref}} = \frac{m_{\text{ring}} \cdot n_{\text{ring}} + m_{\text{lin}} \cdot n_{\text{lin}}}{L_0^3}$$

$$L_{\text{stage}} = \left(\frac{M_{\text{total}}}{\phi_{\text{ref}}}\right)^{1/3}$$

### Linking model

For each existing ring $t$ with length $n_{\text{ring},t}$, the per-target linking intensity upon a cyclisation of length $l_{\text{cyc}}$ is:

$$\mu = A \cdot \frac{n_{\text{ring},t} \cdot l_{\text{cyc}}}{V_{\text{box}}}$$

A link is added with probability $p = 1 - e^{-\mu}$ (Bernoulli trial per target).

### Rate constants

The DSMC engine uses Smoluchowski coagulation ($k_1$) and cyclisation ($k_2$) rates with ideal-chain scaling ($\nu = 0.5$). The ratio $\kappa = k_2/k_1$ controls the cyclised-length distribution; absolute rates only set timescales. Defaults are fitted via JS-divergence minimisation against the MD cyclised-length PMF at `nlin=64`:

| Parameter | Default | Source |
|---|---|---|
| `k1` | 1.0 | Reference scale |
| `k2` | 12840.85 | PMF fit (`fit_to_md.py`) |
| `A` | 0.2093 | Valence fit (`fit_valence_model.py`) |

---

## Setup

```bash
conda create -n newbase python=3.10 -y
conda activate newbase
pip install numpy pandas scipy matplotlib networkx tqdm
```

---

## Repository structure

```
LigationMC/
├── simulation/                    # Core simulation engine and fitting
│   ├── __init__.py
│   ├── main_ligmc.py              # CLI entry point for staged simulations
│   ├── dsmc_engine.py             # Gillespie SSA: merge + cyclisation events
│   ├── network_builder.py         # Ring graph growth and Poisson linking
│   ├── polymer_utils.py           # Physics constants, c* scaling, helpers
│   ├── analysis.py                # Result serialization and plot helpers
│   ├── fit_to_md.py               # Fit k1/k2 ratio from cyclised-length PMF
│   ├── fit_valence_model.py       # Fit linking prefactor A from MD valence data
│   ├── rate_fitting.py            # Fit k1, k2 from MD average-length trajectories
│   └── sweep_gelation.py          # Parameter sweep over (nring, nlin) grid
│
├── plotting/                      # Post-processing visualisation
│   ├── __init__.py
│   ├── compare_gel_point_time.py   # Scatter: sim vs MD gel-point stages
│   ├── compare_links_per_stage.py  # Observed vs predicted links per event
│   ├── plot_sim_vs_md_by_nlin.py   # PMF overlay plots per nlin
│   └── plot_gelation_phase_diagram.py  # 2×2 gelation phase diagram
│
├── DSMC/                          # Legacy MATLAB/Python fitting references
├── run_systems.sh                 # Batch runner for the 16-system MD matrix
├── fit.md                         # Fitting methodology documentation
├── PLAN.md                        # Original project context
└── README.md
```

---

## Usage

### Run a single system

```bash
python -m simulation.main_ligmc \
    --mring 78 --nring 512 --mlin 6 --nlin 64 \
    --n_stages 100 --trials 1000 --workers 32 \
    --nu 0.5 --L 80
```

Key CLI arguments:

| Argument | Description | Default |
|---|---|---|
| `--mring` | Number of initial ring polymers | — |
| `--nring` | Ring polymer length (monomers) | — |
| `--mlin` | Number of linears injected per stage | — |
| `--nlin` | Linear polymer length (monomers) | — |
| `--n_stages` | Number of growth stages | 100 |
| `--trials` | Independent MC trials | 1000 |
| `--workers` | Parallel workers | 1 |
| `--L` | Reference box side length (σ) | 80 |
| `--k1 / --k2` | DSMC rate constants | 1.0 / 12840.85 |
| `--val_A` | Linking prefactor | 0.2093 |
| `--nu` | Flory exponent | 0.5 |
| `--skip_plots` | Skip plot generation | False |

### Run the full 16-system MD matrix

```bash
bash run_systems.sh
```

### Sweep over (nring, nlin) grid

Compute `(mring, mlin)` from concentration constraints (rings at 5c\*, linears at 0.05c\*) and run all systems:

```bash
python -m simulation.sweep_gelation \
    --nring_min 256 --nring_max 1040 --nring_step 16 \
    --nlin_min 64  --nlin_max 224  --nlin_step 16 \
    --ccsr 5.0 --ccsl 0.05 --L 80 \
    --trials 1000 --n_stages 100 --workers 32 \
    --skip_plots \
    --out_dir /path/to/results
```

This produces `sweep_summary.csv` with one row per system.

### Plot the gelation phase diagram

```bash
python -m plotting.plot_gelation_phase_diagram \
    --sweep_csv /path/to/sweep_summary.csv \
    --out_dir /path/to/results
```

Produces a 2×2 panel figure:

| Panel | X | Y | Color | Overlays |
|---|---|---|---|---|
| Top-left | λ × mlin | Stages to 50% | nlin | Per-nlin power-law fits |
| Top-right | nlin | Power-law slope | — | Global power-law fit |
| Bottom-left | nlin | Stages to 50% | nring | Scatter |
| Bottom-right | nring | Stages to 50% | nlin | Per-nlin power-law fits |

---

## Fitting workflows

### 1. Cyclisation PMF fit (k₁, k₂)

Fits the ratio κ = k₂/k₁ by matching the Gillespie SSA cyclised-length PMF to the MD PMF at `nlin=64`. Uses Jensen-Shannon divergence as the objective, with `k1=1` fixed and `k2` optimised in log-space via Nelder-Mead.

```bash
python -m simulation.fit_to_md --nlin 64 --mlin 6 --nu 0.5 --n_trials 5000
```

**Target PMF** (from MD, nlin=64):

| Length | 64 | 128 | 192 | 256 | 320 | 384 |
|---|---|---|---|---|---|---|
| PMF | 0.780 | 0.156 | 0.037 | 0.022 | 0.005 | 0.002 |

**Outputs**: `fitted_k1_k2.json`, `fit_to_md_pmf_comparison.png`

### 2. Linking valence fit (A)

Fits the Poisson linking prefactor A from MD valence data (`summary_all_systems_links_by_size.csv`) using weighted least squares with B=1 fixed.

Model: λ = A · (nring · mring · l_cyc) / V_box

```bash
python -m simulation.fit_valence_model \
    --csv /path/to/summary_all_systems_links_by_size.csv \
    --max_linear_size 256
```

**Outputs**: `fitted_valence_model.json`, `fit_valence_overlay_by_nring.png`

### 3. Rate fitting from MD trajectories (alternative)

Fits k₁ and k₂ by matching average-length time series from MD replicas to a Smoluchowski forward model:

```bash
python -m simulation.rate_fitting \
    --data_root /path/to/replicas \
    --n_molecules 200 --monomer_length 174 --volume 41781609.5
```

**Outputs**: `fitted_rates.json`

### 4. Validation (not fitting)

`plot_sim_vs_md_by_nlin.py` reads fitted k₁/k₂ and generates PMF overlays for multiple nlin values — it does **not** re-fit any parameters.

```bash
python -m plotting.plot_sim_vs_md_by_nlin --nlins 64,96,128,160 --n_trials 5000
```

---

## Data flow

```
┌─────────────────────────┐
│  MD simulation data     │
│  (topology histories)   │
└──────────┬──────────────┘
           │
    ┌──────▼──────┐       ┌───────────────┐
    │ fit_to_md   │       │ rate_fitting   │
    │ (PMF → k₂)  │       │ (ODE → k₁,k₂) │
    └──────┬──────┘       └───────┬───────┘
           │                      │
    ┌──────▼──────────────────────▼──┐
    │       fitted_k1_k2.json        │
    └──────────────┬─────────────────┘
                   │
    ┌──────────────▼──────┐     ┌──────────────────┐
    │ fit_valence_model   │     │ MD valence CSV    │
    │ (WLS → A)           │◄────┤                  │
    └──────┬──────────────┘     └──────────────────┘
           │
    ┌──────▼──────────────────────┐
    │ main_ligmc / sweep_gelation │
    │ (staged simulation)        │
    └──────┬──────────────────────┘
           │
    ┌──────▼──────────────────────┐
    │ plotting scripts            │
    │ (phase diagrams, overlays)  │
    └─────────────────────────────┘
```

---

## Outputs

By default, results are written to:

```
/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc
```

Generated artefacts per system:
- `results_L{L}_mring{M}_nring{N}_mlin{m}_nlin{n}.pkl` — per-trial results
- `results_all.pkl` — consolidated across systems
- `event_timeline.pkl` — full event log
- `gelation_curves.png` — largest-component fraction vs stage
- `degree_distributions.png` — final degree histograms
- `stages_vs_nlin_nring.png` — gelation stage scatter plots

---

## Debugging tips

- Confirm `event.ring_length` is set and non-zero (this is the cyclised linear length).
- If nothing links, check the scale of `A` and `box_volume` — `p` can be near 0 or 1.
- Use deterministic RNG with `NetworkBuilder(seed=...)` for reproducible runs.
- Check `box_length_per_stage` in results to verify the box is rescaling correctly.
- Quick debug print inside `NetworkBuilder.process_cyclisation`:
  ```python
  print('nring_target', nring_target, 'nlin', l_cyc, 'box_volume', box_volume)
  print('mu', mu, 'p', 1 - math.exp(-mu))
  ```
