#!/bin/bash

python ../../src/plumed/main.py \
    --filepath "../../data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_alpha.pdb" \
    --device_index "nan" \
    --device "cuda" \
    --timestep 2 \
    --mdtime 1 \
    --output_dir "../../data/241010_FoldingUponBinding/output/241029/A-synuclein/debug"
