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
    local metad_pace=$3
    local output_dir=$4
    local restart=$5
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,METAD_PACE=$metad_pace,OUTPUT_DIR=$output_dir,RESTART=$restart" $TEMPLATE
}

RESTART_DATE="241204-MetaD2"
DATE="241204-MetaD3"
SYSTEM="Z1-B50L10W"
PROJECT_DIR="/home/mmm1486/projects/openmm-md/data/241109_INFconstruct/output/"


# cmap, d are the cvs in that order
# submit_simulation "${PROJECT_DIR}/${SYSTEM}/241122-Explore/${SYSTEM}_equilibrated.pdb" \
#                 "${SYSTEM}" \
#                 50 \
#                 "${PROJECT_DIR}/data/241109_INFconstruct/output/${SYSTEM}/${DATE}-MetaDPace50"

# submit_simulation "${PROJECT_DIR}/${SYSTEM}/241122-Explore/${SYSTEM}_equilibrated.pdb" \
#                 "${SYSTEM}" \
#                 100 \
#                 "${PROJECT_DIR}/data/241109_INFconstruct/output/${SYSTEM}/${DATE}-MetaDPace100"

submit_simulation "/home/mmm1486/projects/openmm-md/data/241109_INFconstruct/output/${SYSTEM}/${RESTART_DATE}/${SYSTEM}_equilibrated.cif" \
                "${SYSTEM}" \
                500 \
                "/home/mmm1486/projects/openmm-md/data/241109_INFconstruct/output/${SYSTEM}/${DATE}" \
                True


# submit_simulation "${PROJECT_DIR}/${SYSTEM}/241122-Explore/${SYSTEM}_equilibrated.pdb" \
#                 "${SYSTEM}" \
#                 56 \
#                 0.24 \
#                 0.05 \
#                 0 \
#                 103 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/ASYN-A/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/ASYN-G/241122-Explore/A-synuclein_general_equilibrated.pdb" \
#                 "ASYN-G" \
#                 12 \
#                 0.31 \
#                 0.25 \
#                 0 \
#                 88 \
#                 0 \
#                 5 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/ASYN-G/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-A/241122-Explore/CD28_alpha_equilibrated.pdb" \
#                 "CD28-A" \
#                 40 \
#                 0.18 \
#                 0.02 \
#                 0 \
#                 45 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-A/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-B/241122-Explore/CD28_beta_equilibrated.pdb" \
#                 "CD28-B" \
#                 40 \
#                 0.23 \
#                 0.01 \
#                 0 \
#                 41 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-B/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/241122-Explore/CD28_general_equilibrated.pdb" \
#                 "CD28-G" \
#                 48 \
#                 0.11 \
#                 0.27 \
#                 0 \
#                 31 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-G/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-P/241122-Explore/CD28_partial_equilibrated.pdb" \
#                 "CD28-P" \
#                 48 \
#                 0.18 \
#                 0.02 \
#                 0 \
#                 43 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/CD28-P/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/P53-1/241122-Explore/p53_1_equilibrated.pdb" \
#                 "P53-1" \
#                 40 \
#                 0.12 \
#                 0.01 \
#                 0 \
#                 33 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/P53-1/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/P53-2/241122-Explore/p53_2_equilibrated.pdb" \
#                 "P53-2" \
#                 40 \
#                 0.17 \
#                 0.01 \
#                 0 \
#                 59 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/P53-2/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/P53-E/241122-Explore/p53_end_equilibrated.pdb" \
#                 "P53-E" \
#                 40 \
#                 0.27 \
#                 0.05 \
#                 0 \
#                 112 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/P53-E/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/SUMO-1A/241122-Explore/sumo1_equilibrated.pdb" \
#                 "SUMO-1A" \
#                 52 \
#                 0.18 \
#                 0.06 \
#                 0 \
#                 63 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/SUMO-1A/${DATE}-MetaD"

# submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/SUMO-1C/241122-Explore/sumo1c_equilibrated.pdb" \
#                 "SUMO-1C" \
#                 32 \
#                 0.23 \
#                 0.01 \
#                 0 \
#                 63 \
#                 0 \
#                 7 \
#                 "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/SUMO-1C/${DATE}-MetaD"
