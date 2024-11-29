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
    local biasfactor=$3
    local sigma_cv1=$4
    local sigma_cv2=$5
    local grid_min_cv1=$6
    local grid_max_cv1=$7
    local grid_min_cv2=$8
    local grid_max_cv2=$9
    local output_dir=${10}
    echo "Submitting $system with template $TEMPLATE"

    qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,BIASFACTOR=$biasfactor,SIGMA_CV1=$sigma_cv1,SIGMA_CV2=$sigma_cv2,GRID_MIN_CV1=$grid_min_cv1,GRID_MAX_CV1=$grid_max_cv1,GRID_MIN_CV2=$grid_min_cv2,GRID_MAX_CV2=$grid_max_cv2,OUTPUT_DIR=$output_dir" $TEMPLATE
}

DATE="241128"
PROJECT_DIR="/home/mmm1486/projects/openmm-md"


# cmap, d are the cvs in that order
submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/ASYN-A/241122-Explore/A-synuclein_alpha_equilibrated.pdb" \
                "ASYN-A" \
                56 \
                0.24 \
                0.05 \
                0 \
                103 \
                0 \
                7 \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/ASYN-A/${DATE}-MetaD"

submit_simulation "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/ASYN-G/241122-Explore/A-synuclein_general_equilibrated.pdb" \
                "ASYN-G" \
                12 \
                0.31 \
                0.25 \
                0 \
                88 \
                0 \
                5 \
                "${PROJECT_DIR}/data/241010_FoldingUponBinding/output/ASYN-G/${DATE}-MetaD"

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
