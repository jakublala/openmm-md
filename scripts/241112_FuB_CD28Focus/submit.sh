#!/bin/bash

submit_simulation() {
    local input_path=$1
    local output_dir=$2
    local job_name=$3
    local timestep=${4:-2}
    local mdtime=${5:-500}
    local restart_rfile=${6:-None}
    local barrier=${7:-200}
    local padding=${8:-5}
    local upper_wall=${9:-5}

    qsub -N "$job_name" \
        -v "INPUT_PATH=$input_path,OUTPUT_DIR=$output_dir,TIMESTEP=$timestep,MDTIME=$mdtime,RESTART_RFILE=$restart_rfile,BARRIER=$barrier,PADDING=$padding,UPPER_WALL=$upper_wall" \
        template.pbs
}

# CD28-G-UW5P2 (5nm upper wall, 2nm padding)
submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                 "data/241010_FoldingUponBinding/output/CD28-G-UW5P2/241115" \
                 "cd28_uw5p2" \
                 2 \
                 500 \
                 "None" \
                 200 \
                 2 \
                 5

# CD28-G-UW5P4 (5nm upper wall, 4nm padding)
submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                 "data/241010_FoldingUponBinding/output/CD28-G-UW5P4/241115" \
                 "cd28_uw5p4" \
                 2 \
                 500 \
                 "None" \
                 200 \
                 4 \
                 5

# CD28-G-UW8P4 (8nm upper wall, 4nm padding)
submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                 "data/241010_FoldingUponBinding/output/CD28-G-UW8P4/241115" \
                 "cd28_uw8p4" \
                 2 \
                 500 \
                 "None" \
                 200 \
                 4 \
                 8

# python ../../src/plumed/main.py \
#     --filepath "../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
#     --device_index "nan" \
#     --device "cuda" \
#     --timestep 2 \
#     --mdtime 500 \
#     --output_dir "../../data/241010_FoldingUponBinding/output/241112/CD28/alpha_1" \
#     --barrier 200 \
#     --padding 5 \
#     --upper_wall_at 5