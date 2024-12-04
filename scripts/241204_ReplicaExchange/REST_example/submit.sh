#!/bin/bash

#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m n
#$ -l cuda
#$ -l gpu=1
#$ -pe one_node 1
#$ -l gtx1080ti
#$ -P dev_GPU


export PATH=/nfs/working/deep_learn/xiaowei/miniconda3/bin:$PATH

python REST2_example.py > REST2_example.out
