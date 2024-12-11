#!/bin/bash

# Get the template PBS file and set PROJECT_DIR
if [[ "$1" == "-f" && -n "$2" ]]; then
    TEMPLATE="$2"
    if [[ ! -f "$TEMPLATE" ]]; then
        echo "Specified template file $TEMPLATE not found"
        exit 1
    fi
    # Set PROJECT_DIR based on template name
    if [[ "$TEMPLATE" == "hx1.pbs" ]]; then
        PROJECT_DIR="/gpfs/home/jl24018/projects/openmm-md"
    elif [[ "$TEMPLATE" == "mmm.pbs" ]]; then
        PROJECT_DIR="/home/mmm1486/projects/openmm-md"
    fi
else
    if [[ -f "hx1.pbs" ]]; then
        PROJECT_DIR="/gpfs/home/jl24018/projects/openmm-md"
        TEMPLATE="hx1.pbs"
    elif [[ -f "cx3.pbs" ]]; then
        TEMPLATE="cx3.pbs" 
    elif [[ -f "mmm.pbs" ]]; then
        PROJECT_DIR="/home/mmm1486/projects/openmm-md"
        TEMPLATE="mmm.pbs"
    else
        echo "No PBS template file found"
        exit 1
    fi
fi


submit_simulation() {
    local filepath=$1
    local system=$2
    local output_dir=$3
    local device_precision=$4
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,OUTPUT_DIR=$output_dir,DEVICE_PRECISION=$device_precision" $TEMPLATE
}

echo ${PROJECT_DIR}

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                "CD28-G-double" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241209-NoPBC-double" \
                "double"


submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                "CD28-G-mixed" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241210-NoPBC-mixed" \
                "mixed"


submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                "CD28-G-single" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241210-NoPBC-single" \
                "single"


