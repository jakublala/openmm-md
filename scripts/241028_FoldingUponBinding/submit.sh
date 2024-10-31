#!/bin/bash

submit_simulation() {
    local input_path=$1
    local output_dir=$2
    local job_name=$3
    local timestep=${4:-2}
    local mdtime=${5:-500}

    qsub -N $job_name -v INPUT_PATH="$input_path",OUTPUT_DIR="$output_dir",TIMESTEP="$timestep",MDTIME="$mdtime" template.pbs
}

# Example usage:
# submit_simulation "data/241010_FoldingUponBinding/input/p53/p53_1.pdb" \
#                  "data/241010_FoldingUponBinding/output/${date}/p53" \
#                  "p53_sim" \
#                  2 \
#                  1

date="241029"


# A-synuclein simulations
for i in {1..5}; do
    submit_simulation "data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_alpha.pdb" \
                     "data/241010_FoldingUponBinding/output/${date}/A-synuclein/alpha_$i" \
                     "asyn_alpha_$i" \
                     2 \
                     500

    # submit_simulation "data/241010_FoldingUponBinding/input/A-synuclein/A-synuclein_general.pdb" \
    #                  "data/241010_FoldingUponBinding/output/${date}/A-synuclein/general_$i" \
    #                  "asyn_general_$i" \
    #                  2 \
    #                  500
done

# # CD28 simulations
# for i in {1..5}; do
#     submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_alpha.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/CD28/alpha_$i" \
#                      "cd28_alpha_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_beta.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/CD28/beta_$i" \
#                      "cd28_beta_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/CD28/general_$i" \
#                      "cd28_general_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/CD28/CD28_partial.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/CD28/partial_$i" \
#                      "cd28_partial_$i" \
#                      2 \
#                      500
# done

# # p53 simulations
# for i in {1..5}; do
#     submit_simulation "data/241010_FoldingUponBinding/input/p53/p53_1.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/p53/1_$i" \
#                      "p53_1_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/p53/p53_2.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/p53/2_$i" \
#                      "p53_2_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/p53/p53_end.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/p53/end_$i" \
#                      "p53_end_$i" \
#                      2 \
#                      500
# done

# # SUMO simulations
# for i in {1..5}; do
#     submit_simulation "data/241010_FoldingUponBinding/input/SUMO/sumo1.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/SUMO/1_$i" \
#                      "sumo1_$i" \
#                      2 \
#                      500

#     submit_simulation "data/241010_FoldingUponBinding/input/SUMO/sumo1c.pdb" \
#                      "data/241010_FoldingUponBinding/output/${date}/SUMO/1c_$i" \
#                      "sumo1c_$i" \
#                      2 \
#                      500
# done