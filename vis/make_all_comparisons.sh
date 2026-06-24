#!/usr/bin/env bash
# Generate all three final comparison figures (PMF, Valence, Gel Point)

set -e

RESULTS_DIR="./results"
OUT_DIR="${RESULTS_DIR}/final_figures"
mkdir -p "${OUT_DIR}"

echo "1/3: Generating Polymerisation PMF Comparison Grid..."
python3 vis/plot_sim_vs_md_by_nlin.py \
    --nlins 64,96,128,160 \
    --fit_json parameters/fitted_k1_k2.json \
    --out_dir "${OUT_DIR}" \
    --n_trials 5000

echo ""
echo "2/3: Generating Valence Linking Parity Plot..."
python3 vis/plot_valence_md_comparison.py \
    --md-csv /storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels/results/histories/summary_all_systems_links_by_size_s0.csv \
    --model-json parameters/fitted_valence_model.json \
    --sim-root /storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc/md_runs_dsmc \
    --out-png "${OUT_DIR}/valence_comparison.png"

echo ""
echo "3/3: Generating Gel Point Parity Plot..."
python3 vis/compare_gel_point_time.py \
    --gel-file /storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels/values/network/gel_point_time_avgstd.txt \
    --results-root /storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc/md_runs_dsmc \
    --out-csv "${OUT_DIR}/gel_time_compare.csv" \
    --out-png "${OUT_DIR}/gel_time_compare.png"

echo ""
echo "Done! All final figures generated in: ${OUT_DIR}"
ls -la "${OUT_DIR}"
