#!/bin/bash

submit_simulation() {
    local input_path=$1
    local output_dir=$2
    local job_name=$3
    local timestep=${4:-2}
    local mdtime=${5:-500}
    local restart_rfile=${6:-None}

    qsub -N $job_name -v INPUT_PATH="$input_path",OUTPUT_DIR="$output_dir",TIMESTEP="$timestep",MDTIME="$mdtime",RESTART_RFILE="$restart_rfile" template.pbs
}

# A-synuclein simulations
for i in {1..5}; do
    submit_simulation "data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_alpha.pdb" \
                     "data/241010_FoldingUponBinding/output/241028/A-synuclein/alpha_$i" \
                     "asyn_alpha_$i" \
                     2 \
                     500 \
                     "../../data/241010_FoldingUponBinding/output/241029/A-synuclein/alpha_$i/A-synuclein_alpha.state"

    submit_simulation "data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_general.pdb" \
                     "data/241010_FoldingUponBinding/output/241028/A-synuclein/general_$i" \
                     "asyn_general_$i" \
                     2 \
                     500 \
                     "../../data/241010_FoldingUponBinding/output/241029/A-synuclein/general_$i/A-synuclein_general.state"
done

# CD28 simulations
for i in {1..5}; do
    submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_alpha.pdb" \
                     "data/241010_FoldingUponBinding/output/241028/CD28/alpha_$i" \
                     "cd28_alpha_$i" \
                     2 \
                     500 \
                     "../../data/241010_FoldingUponBinding/output/241029/CD28/alpha_$i/CD28_alpha.state"

    submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_beta.pdb" \
                     "data/241010_FoldingUponBinding/output/241028/CD28/beta_$i" \
                     "cd28_beta_$i" \
                     2 \
                     500 \
                     "../../data/241010_FoldingUponBinding/output/241029/CD28/beta_$i/CD28_beta.state"

    submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
                     "data/241010_FoldingUponBinding/output/241028/CD28/general_$i" \
                     "cd28_general_$i" \
                     2 \
                     500 \
                     "../../data/241010_FoldingUponBinding/output/241029/CD28/general_$i/CD28_general.state"

    submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_partial.pdb" \
                     "data/241010_FoldingUponBinding/output/241028/CD28/partial_$i" \
                     "cd28_partial_$i" \
                     2 \
                     500 \
                     "../../data/241010_FoldingUponBinding/output/241029/CD28/partial_$i/CD28_partial.state"
done

# # p53 simulations
# for i in {1..5}; do
#     submit_simulation "data/241010_FoldingUponBinding/input/p53/p53_1.pdb" \
#                      "data/241010_FoldingUponBinding/output/241028/p53/1_$i" \
#                      "p53_1_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/p53/p53_2.pdb" \
#                      "data/241010_FoldingUponBinding/output/241028/p53/2_$i" \
#                      "p53_2_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/p53/p53_end.pdb" \
#                      "data/241010_FoldingUponBinding/output/241028/p53/end_$i" \
#                      "p53_end_$i" \
#                      2 \
#                      500
# done

# # SUMO simulations
# for i in {1..5}; do
#     submit_simulation "data/241010_FoldingUponBinding/input/SUMO/sumo1.pdb" \
#                      "data/241010_FoldingUponBinding/output/241028/SUMO/1_$i" \
#                      "sumo1_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/SUMO/sumo1c.pdb" \
#                      "data/241010_FoldingUponBinding/output/241028/SUMO/1c_$i" \
#                      "sumo1c_$i" \
#                      2 \
#                      500
# done