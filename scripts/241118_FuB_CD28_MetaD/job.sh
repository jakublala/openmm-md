#!/bin/bash

docker run \
    --gpus all \
    -v $HOME/phd/openmm-md:/app \
    -w /app \
    --user $(id -u):$(id -g) \
    jakublala/openmm-md \
    bash -c "pip install -e . && cd scripts/241118_FuB_CD28_MetaD && python run.py"