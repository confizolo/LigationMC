#!/bin/bash
#

conda activate newbase

systems=(
    # "80 27 1024 6 64 1 50",
    "80 27 1024 3 96 1 50",
    "80 27 1024 2 128 1 50",
    "80 27 1024 1 160 1 50",
    "80 42 768 6 64 1 50",
    "80 42 768 3 96 1 50",
    "80 42 768 2 128 1 50",
    "80 42 768 1 160 1 50",
    "80 78 512 6 64 1 50",
    "80 78 512 3 96 1 50",
    "80 78 512 2 128 1 50",
    "80 78 512 1 160 1 50",
    "80 222 256 6 64 1 50",
    "80 222 256 3 96 1 50",
    "80 222 256 2 128 1 50",
    "80 222 256 1 160 1 50",
)

for system in "${systems[@]}"; do
    echo "Running system configuration: ${system}"
    read -r L mr nr ml nl rep time <<< "$system"

    python main_ligmc.py --progress --workers 30 --L $L --mring $mr --nring $nr --mlin $ml --nlin $nl

done