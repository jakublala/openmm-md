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
# Note: this is a bit problematic, as changes in environment.yml are not reflected alone
# in the image change, so need to manually update the image
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Test installation of OpenMM
RUN python -m openmm.testInstallation

# Test installation of OpenMM and CUDA
RUN python -c "import openmm as mm; print('Available platforms:', [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])"


# Install PLUMED with OPES module
RUN git clone https://github.com/plumed/plumed2 \
    && cd plumed2 \
    && ./configure --prefix=$CONDA_DIR --enable-modules=opes,sasa \
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

# For development, we'll use volume mounting
# Do NOT install the package during build
WORKDIR /app

# Note: When running the container, mount your code directory like:
# docker run -v /path/to/your/code:/app ...