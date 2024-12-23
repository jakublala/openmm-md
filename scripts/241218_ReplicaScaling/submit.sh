#!/bin/bash

# Detect environment based on current working directory
CURRENT_PATH=$(pwd)

if [[ "$CURRENT_PATH" == *"gpfs"* ]]; then
    PROJECT_DIR="/gpfs/home/jl24018/projects/openmm-md"
    TEMPLATE="hx1.pbs"
elif [[ "$CURRENT_PATH" == *"mmm1486"* ]]; then
    # PROJECT_DIR="/home/mmm1486/projects/openmm-md"


    # having some issues with MPI
    PROJECT_DIR="/home/mmm1486/Scratch/openmm-md"

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
    local n_replicas=$4
    local gpu=$5
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -l gpu=$gpu -v "FILEPATH=$filepath,SYSTEM=$system,OUTPUT_DIR=$output_dir,N_REPLICAS=$n_replicas" $TEMPLATE
}

echo ${PROJECT_DIR}

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
                "ReplicaScaling2-GPU1" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling2-GPU1" \
                2 \
                1

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling4-GPU1" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling4-GPU1" \
#                 4 \
#                 1

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling8-GPU1" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling8-GPU1" \
#                 8 \
#                 1



# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling2-GPU2" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling2-GPU2" \
#                 2 \
#                 2

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling4-GPU2" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling8-GPU2" \
#                 4 \
#                 2

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling8-GPU2" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling8-GPU2" \
#                 8 \
#                 2



# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling2-GPU4" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling2-GPU4" \
#                 2 \
#                 4

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling2-GPU4" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling4-GPU4" \
#                 4 \
#                 4

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241211-NewUW/CD28_general_equilibrated.cif" \
#                 "ReplicaScaling8-GPU4" \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241218-ReplicaScaling8-GPU4" \
#                 8 \
#                 4



