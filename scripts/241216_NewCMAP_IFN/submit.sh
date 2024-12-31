#!/bin/bash

# Detect environment based on current working directory
CURRENT_PATH=$(pwd)

if [[ "$CURRENT_PATH" == *"gpfs"* ]]; then
    PROJECT_DIR="/gpfs/home/jl24018/projects/openmm-md"
    TEMPLATE="hx1.pbs"
elif [[ "$CURRENT_PATH" == *"mmm1486"* ]]; then
    PROJECT_DIR="/home/mmm1486/projects/openmm-md"
    TEMPLATE="mmm.pbs"
else
    echo "Unknown environment: $CURRENT_PATH"
    exit 1
fi

echo "Detected environment: Using $TEMPLATE with project dir: ${PROJECT_DIR}"


submit_simulation() {
    local filepath=$1
    local system=$2
    local output_dir=$3
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,OUTPUT_DIR=$output_dir" $TEMPLATE
}

echo ${PROJECT_DIR}
submit_simulation "${PROJECT_DIR}/data/241109_INFconstruct/output/Z1-B50L10W/241216-NewCMAP/Z1-B50L10W_equilibrated.cif" \
                "Z1-B50L10W-CMAP" \
                "${PROJECT_DIR}/data/241109_INFconstruct/output/Z1-B50L10W/241231-LongCMAP" \
