# Use NVIDIA CUDA base image with Ubuntu 22.04
FROM mambaorg/micromamba:git-c7619c9-cuda11.8.0-ubuntu22.04

# Set environment variables
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
ENV LD_LIBRARY_PATH=$CONDA_DIR/lib:$LD_LIBRARY_PATH

# Create environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Enable conda environment during build
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install conda packages
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Install PLUMED with OPES module
RUN git clone https://github.com/plumed/plumed2 \
    && cd plumed2 \
    && ./configure --prefix=$CONDA_DIR --enable-modules=opes \
    && make -j4 \
    && make install \
    && cd .. \
    && rm -rf plumed2

# Install OpenMM-PLUMED plugin
RUN git clone https://github.com/openmm/openmm-plumed \
    && cd openmm-plumed \
    && mkdir build && cd build \
    && cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=$CONDA_DIR \
        -DOPENMM_DIR=$CONDA_DIR \
        -DPLUMED_INCLUDE_DIR=$CONDA_DIR/include/plumed \
        -DPLUMED_LIBRARY_DIR=$CONDA_DIR/lib \
        -DPYTHON_EXECUTABLE=$CONDA_DIR/bin/python \
        -DPIP_EXECUTABLE=$CONDA_DIR/bin/pip \
        -DSWIG_EXECUTABLE=$CONDA_DIR/bin/swig \
        -DPLUMED_BUILD_PYTHON_WRAPPERS=ON \
    && sed -i '1770i#include <memory>' $CONDA_DIR/include/plumed/wrapper/Plumed.h \
    && make -j4 \
    && make install \
    && make PythonInstall \
    && cd ../.. \
    && rm -rf openmm-plumed

# For your package installation, you have two options:
# 1. Install during build (uncomment and modify the following):
# COPY . /app
# RUN pip install -e /app

# 2. Use a volume mount during runtime (recommended for development)
# Then you would mount your code directory when running the container

WORKDIR /app