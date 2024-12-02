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


submit_simulation "../../data/241109_INFconstruct/input/Q7-B30L4.pdb" \
                "Q7-B30L4"

submit_simulation "../../data/241109_INFconstruct/input/Q7-B30L7.pdb" \
                "Q7-B30L7"

submit_simulation "../../data/241109_INFconstruct/input/Q7-B30L10.pdb" \
                "Q7-B30L10"

submit_simulation "../../data/241109_INFconstruct/input/Q7-B40L10.pdb" \
                "Q7-B40L10"

submit_simulation "../../data/241109_INFconstruct/input/Q7-B40L10W.pdb" \
                "Q7-B40L10W"

submit_simulation "../../data/241109_INFconstruct/input/Q7-B50L10W.pdb" \
                "Q7-B50L10W"

submit_simulation "../../data/241109_INFconstruct/input/PQ19-B30L4.pdb" \
                "PQ19-B30L4"

submit_simulation "../../data/241109_INFconstruct/input/PQ19-B30L7.pdb" \
                "PQ19-B30L7"

submit_simulation "../../data/241109_INFconstruct/input/PQ19-B30L10.pdb" \
                "PQ19-B30L10"

submit_simulation "../../data/241109_INFconstruct/input/PQ19-B40L10.pdb" \
                "PQ19-B40L10"

submit_simulation "../../data/241109_INFconstruct/input/PQ19-B40L10W.pdb" \
                "PQ19-B40L10W"

submit_simulation "../../data/241109_INFconstruct/input/PQ19-B50L10W.pdb" \
                "PQ19-B50L10W"

submit_simulation "../../data/241109_INFconstruct/input/Z1-B30L4.pdb" \
                "Z1-B30L4"

submit_simulation "../../data/241109_INFconstruct/input/Z1-B30L7.pdb" \
                "Z1-B30L7"

submit_simulation "../../data/241109_INFconstruct/input/Z1-B30L10.pdb" \
                "Z1-B30L10"

submit_simulation "../../data/241109_INFconstruct/input/Z1-B40L10.pdb" \
                "Z1-B40L10"

submit_simulation "../../data/241109_INFconstruct/input/Z1-B40L10W.pdb" \
                "Z1-B40L10W"

submit_simulation "../../data/241109_INFconstruct/input/Z1-B50L10W.pdb" \
                "Z1-B50L10W"