#!/bin/bash

#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m n
#$ -l cuda
#$ -l gpu=4
#$ -pe one_node 1
#$ -l gtx1080ti
#$ -P dev_GPU


export PATH=/nfs/working/deep_learn/xiaowei/miniconda3/bin:$PATH

export SLURM_JOB_NODELIST=$HOSTNAME
#export CUDA_VISIBLE_DEVICES=0,1
export SLURMD_NODENAME=$HOSTNAME

python build_mpirun_configfile.py --configfilepath configfile --hostfilepath hostfile "python REST2_example.py"
mpiexec.hydra -f hostfile -configfile configfile > REST2_example_mpi.out 
