#!/bin/bash

# Array of all systems
systems=(
    "A-synuclein_alpha"
    "A-synuclein_general"
    "CD28_alpha"
    "CD28_beta"
    "CD28_partial"
    "CD28_general"
    "p53_1"
    "p53_2"
    "p53_end"
    "sumo_1"
    "sumo_1c"
)

# Loop through each system and submit a job
for system in "${systems[@]}"; do
    # Pipe directly to sbatch using process substitution
    sbatch <(sed "s/{system}/$system/g" analysis)
done