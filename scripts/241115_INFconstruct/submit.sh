#!/bin/bash

submit_simulation() {
    local input_path=$1
    local output_dir=$2
    local job_name=$3
    local timestep=${4:-2}
    local mdtime=${5:-500}
    local restart_rfile=${6:-None}

    qsub -N $job_name -v INPUT_PATH="$input_path",OUTPUT_DIR="$output_dir",TIMESTEP="$timestep",MDTIME="$mdtime" template_hx1.pbs
}

# submit_simulation "data/241109_INFconstruct/input/Q7-B30L4.pdb" \
#                  "data/241109_INFconstruct/output/Q7-B30L4/241111" \
#                  "Q7-B30L4" \
#                  2 \
#                  100

for i in "Q7" "Z1"; do
    for j in "B30L4" "B30L7" "B30L10" "B40L10" "B40L10W" "B50W"; do
        if [ $i == "Q7" ] && [ $j == "B30L4" ]; then
            continue
        fi
        submit_simulation "data/241109_INFconstruct/input/$i-$j.pdb" \
                         "data/241109_INFconstruct/output/$i-$j/241111" \
                         "$i-$j" \
                         2 \
                         100
    done
done