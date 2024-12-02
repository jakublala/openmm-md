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
    local metad_pace=$3
    local output_dir=$4
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,METAD_PACE=$metad_pace,OUTPUT_DIR=$output_dir" $TEMPLATE
}

DATE="241202"
PROJECT_DIR="/home/mmm1486/projects/openmm-md"

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241128-MetaD/CD28_general_equilibrated.pdb" \
                "CD28-G" \
                50 \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/${DATE}-Pace50"


submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241128-MetaD/CD28_general_equilibrated.pdb" \
                "CD28-G" \
                100 \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/${DATE}-Pace100"


submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241128-MetaD/CD28_general_equilibrated.pdb" \
                "CD28-G" \
                200 \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/${DATE}-Pace200"