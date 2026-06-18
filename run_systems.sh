#!/bin/bash
# Run the standard MD-matched system matrix with the Julia engine.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="${ROOT_DIR}/results_matrix"

mkdir -p "$OUT_DIR"

systems=(
    "80 27 1024 3 96 1000 100"
    "80 27 1024 2 128 1000 100"
    "80 27 1024 1 160 1000 100"
    "80 42 768 6 64 1000 100"
    "80 42 768 3 96 1000 100"
    "80 42 768 2 128 1000 100"
    "80 42 768 1 160 1000 100"
    "80 78 512 6 64 1000 100"
    "80 78 512 3 96 1000 100"
    "80 78 512 2 128 1000 100"
    "80 78 512 1 160 1000 100"
    "80 222 256 6 64 1000 100"
    "80 222 256 3 96 1000 100"
    "80 222 256 2 128 1000 100"
    "80 222 256 1 160 1000 100"
)

for system in "${systems[@]}"; do
    echo "Running system configuration: ${system}"
    read -r L mr nr ml nl trials n_stages <<< "$system"

    julia --project="${ROOT_DIR}/simulation_jl" "${ROOT_DIR}/simulation_jl/RunSingle.jl" \
        --L "$L" --mring "$mr" --nring "$nr" --mlin "$ml" --nlin "$nl" \
        --trials "$trials" --n_stages "$n_stages" --out_dir "$OUT_DIR"

done