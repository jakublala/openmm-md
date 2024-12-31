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

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241222-LongCMAP/CD28_general_equilibrated.cif" \
                "CD28-G-LongCMAP" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241231-LongCMAP"