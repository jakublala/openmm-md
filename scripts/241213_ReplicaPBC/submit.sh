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
    local t_min=$4
    local t_max=$5
    local n_replicas=$6
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,OUTPUT_DIR=$output_dir,T_MIN=$t_min,T_MAX=$t_max,N_REPLICAS=$n_replicas" $TEMPLATE
}

echo ${PROJECT_DIR}

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
                "CD28-G-ReplicaPBC-300-310-8" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241229-ReplicaPBC-300-310-8" \
                "300" \
                "310" \
                "8"

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
                "CD28-G-ReplicaPBC-300-325-8" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241229-ReplicaPBC-300-325-8" \
                "300" \
                "325" \
                "8"


submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
                "CD28-G-ReplicaPBC-300-350-8" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241229-ReplicaPBC-300-350-8" \
                "300" \
                "350" \
                "8"