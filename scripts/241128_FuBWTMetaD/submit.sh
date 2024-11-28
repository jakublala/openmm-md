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
    local filepath=$1
    local system=$2
    local biasfactor=$3
    local sigma=$4
    local grid_min=$5
    local grid_max=$6

    qsub -N "$system" \
        -v "FILEPATH=$filepath,SYSTEM=$system,BIASFACTOR=$biasfactor,SIGMA=$sigma,GRID_MIN=$grid_min,GRID_MAX=$grid_max" \
        $TEMPLATE
}


# cmap, d are the cvs in that order
submit_simulation "../../data/241010_FoldingUponBinding/output/ASYN-A/241122-Explore/A-synuclein_alpha_equilibrated.pdb" \
                "ASYN-A" \
                56 \
                "0.24,0.05" \
                "0,0" \
                "90,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_general.pdb" \
                "ASYN-G" \
                12 \
                "0.31,0.25" \
                "0,0" \
                "60,5"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_alpha.pdb" \
                "CD28-A" \
                40 \
                "0.18,0.02" \
                "0,0" \
                "40,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_beta.pdb" \
                "CD28-B" \
                40 \
                "0.23,0.01" \
                "0,0" \
                "30,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                "CD28-G" \
                48 \
                "0.11,0.27" \
                "0,0" \
                "30,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_partial.pdb" \
                "CD28-P" \
                48 \
                "0.18,0.02" \
                "0,0" \
                "40,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/p53/p53_1.pdb" \
                "P53-1" \
                40 \
                "0.12,0.01" \
                "0,0" \
                "35,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/p53/p53_2.pdb" \
                "P53-2" \
                40 \
                "0.17,0.01" \
                "0,0" \
                "50,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/p53/p53_end.pdb" \
                "P53-E" \
                40 \
                "0.27,0.05" \
                "0,0" \
                "100,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/SUMO/sumo1.pdb" \
                "SUMO-1A" \
                52 \
                "0.18,0.06" \
                "0,0" \
                "60,7"

submit_simulation "../../data/241010_FoldingUponBinding/input/SUMO/sumo1c.pdb" \
                "SUMO-1C" \
                32 \
                "0.23,0.01" \
                "0,0" \
                "60,7"