#!/bin/bash
source activate openmm

# Set OpenMM CPU threads to match PBS allocation
export OPENMM_CPU_THREADS=256

filepath='data/241010_FoldingUponBinding/input/p53/p53_1.pdb'

python ../../src/plumed/main.py \
    --filepath "../../$filepath" \
    --timestep 2 \
    --mdtime 2 \
    --device_index 'nan' \
    --device "cpu" \
    --output_dir "../../data/241010_FoldingUponBinding/output/241114/p53_CPU_workstation" \
    --split_chains True \
    --logging_frequency 1 \
    --plumed_type "opes" \