#!/bin/bash -l

# Request a number of GPU cards, in this case 2 (the maximum)
#$ -l gpu=1

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:10:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=80G

# Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G

# Set the name of the job.
#$ -N GPUJob

# Set the working directory to somewhere in your scratch space.
#$ -wd /home/mmm1486/Scratch/openmm-md

#$ -P Free
#$ -A Imperial_Mat

# Change into temporary directory to run work
cd $TMPDIR

module load apptainer
module load cuda/11.8

apptainer pull --force openmm-md.sif docker://jakublala/openmm-md:latest

apptainer exec openmm-md.sif bash -c "pip install -e $HOME/projects/openmm-md && python $HOME/projects/openmm-md/src/plumed/main.py \
    --filepath 'data/241109_INFconstruct/input/PQ19/PQ19_B30L4.pdb' \
    --device_index 'nan' \
    --device 'cuda' \
    --timestep 2 \
    --mdtime 5 \
    --output_dir 'data/241111_DockerApptainer/output'"

# 10. Preferably, tar-up (archive) all output files onto the shared scratch area
tar zcvf $HOME/Scratch/files_from_job_$JOB_ID.tar.gz $TMPDIR