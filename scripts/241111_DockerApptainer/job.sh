#!/bin/bash -l

# Request a number of GPU cards, in this case 2 (the maximum)
#$ -l gpu=1

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:10:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=1G

# Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=15G

# Set the name of the job.
#$ -N GPUJob

# Set the working directory to somewhere in your scratch space.
#$ -wd /home/mmm1486/Scratch/openmm-md

# Change into temporary directory to run work
cd $TMPDIR

module load apptainer
module load cuda/11.8
apptainer pull docker://jakublala/openmm-md:latest

# Run the application - the line below is just a random example.
mygpucode

# 10. Preferably, tar-up (archive) all output files onto the shared scratch area
tar zcvf $HOME/Scratch/files_from_job_$JOB_ID.tar.gz $TMPDIR


module load anaconda3/personal
source activate openmm
module load cuda/11.4.2

# Read input parameters
INPUT_PATH=${INPUT_PATH}
OUTPUT_DIR=${OUTPUT_DIR}
TIMESTEP=${TIMESTEP:-2}      # Default to 2 if not specified
MDTIME=${MDTIME:-500}        # Default to 500 if not specified
RESTART_RFILE=${RESTART_RFILE:-None}

python ../../src/plumed/main.py \
    --filepath "../../$INPUT_PATH" \
    --device_index "nan" \
    --device "cuda" \
    --timestep $TIMESTEP \
    --mdtime $MDTIME \
    --restart_rfile $RESTART_RFILE \
    --output_dir "../../$OUTPUT_DIR"