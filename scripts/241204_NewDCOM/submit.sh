#!/bin/bash

# Get the template PBS file
if [[ "$1" == "-f" && -n "$2" ]]; then
    TEMPLATE="$2"
    if [[ ! -f "$TEMPLATE" ]]; then
        echo "Specified template file $TEMPLATE not found"
        exit 1
    fi
else
    if [[ -f "hx1.pbs" ]]; then
        TEMPLATE="hx1.pbs"
    elif [[ -f "cx3.pbs" ]]; then
        TEMPLATE="cx3.pbs" 
    elif [[ -f "mmm.pbs" ]]; then
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
    local restart=$4
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,OUTPUT_DIR=$output_dir,RESTART=$restart" $TEMPLATE
}

PROJECT_DIR="/home/mmm1486/projects/openmm-md"

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241204-Long/CD28_general_equilibrated.pdb" \
                "CD28-G" \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241204-Long2" \
                True
