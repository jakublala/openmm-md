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

    qsub -N "$system" \
        -v "FILEPATH=$filepath,SYSTEM=$system" \
        $TEMPLATE
}


submit_simulation "../../data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_alpha.pdb" \
                "ASYN-A"

submit_simulation "../../data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_general.pdb" \
                "ASYN-G"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_alpha.pdb" \
                "CD28-A"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_beta.pdb" \
                "CD28-B"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                "CD28-G"

submit_simulation "../../data/241010_FoldingUponBinding/input/CD28/CD28_partial.pdb" \
                "CD28-P"

submit_simulation "../../data/241010_FoldingUponBinding/input/p53/p53_1.pdb" \
                "P53-1"

submit_simulation "../../data/241010_FoldingUponBinding/input/p53/p53_2.pdb" \
                "P53-2"

submit_simulation "../../data/241010_FoldingUponBinding/input/p53/p53_end.pdb" \
                "P53-E"


submit_simulation "../../data/241010_FoldingUponBinding/input/SUMO/sumo1.pdb" \
                "SUMO-1A"

submit_simulation "../../data/241010_FoldingUponBinding/input/SUMO/sumo1c.pdb" \
                "SUMO-1C"