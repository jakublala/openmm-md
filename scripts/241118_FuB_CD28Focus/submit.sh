#!/bin/bash

# Get the template PBS file
if [[ -f "hx1.pbs" ]]; then
    TEMPLATE="hx1.pbs"
elif [[ -f "cx3.pbs" ]]; then
    TEMPLATE="cx3.pbs" 
else
    echo "No PBS template file found"
    exit 1
fi


submit_simulation() {
    local system=$1
    local padding=${2:-5}
    local upper_wall=${3:-5}

    qsub -N "$system" \
        -v "SYSTEM=$system,PADDING=$padding,UPPER_WALL=$upper_wall" \
        $TEMPLATE
}

# # CD28-G-UW5P2 (5nm upper wall, 2nm padding)
# submit_simulation "CD28-G-UW5P2" \
#                 2 \
#                 5

# # CD28-G-UW5P4 (5nm upper wall, 4nm padding)
# submit_simulation "CD28-G-UW5P4" \
#                 4 \
#                 5

# CD28-G-UW8P4 (8nm upper wall, 4nm padding)
submit_simulation "CD28-G-UW8P4" \
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