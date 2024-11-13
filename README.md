# openmm-md

## Install

## Create environment with OpenMM and packages.
Create env.
```
conda create -n openmm
conda install python pdbfixer
conda install -c conda-forge openmm
pip install -e .
```

Run `pip install -e .`
Then also run `conda install cudatoolkit=11.8` based on your CUDA version to speed up OpenMM.

## Install PLUMED from source
```
git clone https://github.com/plumed/plumed2
cd plumed2
./configure --prefix=/rds/general/user/jl24018/home/anaconda3/envs/openmm --enable-modules=opes
make -j 4
make install
```

## Install OpenMM-Plumed plugin from source

``
git clone https://github.com/openmm/openmm-plumed
git clone https://github.com/openmm/openmm-plumed
cd openmm-plumed
mkdir build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/rds/general/user/jl24018/home/anaconda3/envs/openmm \
  -DOPENMM_DIR=/rds/general/user/jl24018/home/anaconda3/envs/openmm \
  -DPLUMED_INCLUDE_DIR=/rds/general/user/jl24018/home/anaconda3/envs/openmm/include/plumed \
  -DPLUMED_LIBRARY_DIR=/rds/general/user/jl24018/home/anaconda3/envs/openmm/lib \
  -DPYTHON_EXECUTABLE=/rds/general/user/jl24018/home/anaconda3/envs/openmm/bin/python \
  -DPIP_EXECUTABLE=/rds/general/user/jl24018/home/anaconda3/envs/openmm/bin/pip \
  -DSWIG_EXECUTABLE=/rds/general/user/jl24018/home/anaconda3/envs/openmm/bin/swig \
  -DPLUMED_BUILD_PYTHON_WRAPPERS=ON
``
Add #include <memory> to file "conda-folder"/envs/openmm/include/plumed/wrapper/Plumed.h to line 1770, before running
```
make -j4
make install
conda install swig
make PythonInstall
```
Then need to apply to the PBS script:
```
export LD_LIBRARY_PATH=/rds/general/user/jl24018/home/anaconda3/envs/openmm/lib:$LD_LIBRARY_PATH
```



```
plumed
py-plumed
openmm-plumed
```
Need to be installed from source.
Instructions in Notion.


Python wrapper didn't seem to install properly, so had to do:
`pip install plumed` and then set the kernel:
```
import os
os.environ["PLUMED_KERNEL"] = "/home/jakub/anaconda3/envs/openmm/lib/libplumedKernel.so"
import plumed
p=plumed.Plumed()
```

Also need to do this`
```
export LD_LIBRARY_PATH="/home/jakub/anaconda3/envs/openmm/lib:$LD_LIBRARY_PATH"
```


### Required packages
- `openmm cudatoolkit=11.8`
- `fire`
- `mdtraj`
- `pdbfixer`
- `biopython`


- `netCDF4` - speeds up MDanalysis logging of the trajectory - maybe doesn't?? and slows down actually all of it!?! need to check!!!


### Getting data from HX1 to Workstation
`scp -r jl24018@login.hx1.hpc.ic.ac.uk:/gpfs/home/jl24018/projects/openmm-md/scripts/241010_FoldingUponBinding/output/output.zip output.zip`


### Getting data from rclone / onedrive
```
cd openmm-md
rclone sync --ignore-size --retries 5 --low-level-retries 15 onedrive:data/241010_FoldingUponBinding data/241010_FoldingUponBinding
```

### What should be in the metadata config of a simulation?
- `type`: OPES or MD or fixed bias MD
- `temp`: temperature
- `timestep`: time step in fs
- `padding`: how large is the water padding
- `checkpoint frequency`
- `logging frequency`
- `OPES things, PLUMED things`
### Dockerfile
Building docker image.
```
docker build -t jakublala/openmm-md:latest .
docker push jakublala/openmm-md:latest
```
Use it in a jobscript with PBS job scheduler.
Note that the current (latest) image supports CUDA 11.8! I am not sure if that has to be the same as the one loaded? Maybe not... :/
```
module load apptainer
apptainer pull docker://jakublala/openmm-md:latest
```

Run interactively (debugging and development):
```
docker run -it jakublala/openmm-md
```


Check some python code in Apptainer
```
apptainer exec openmm-md_v1.0.0.sif python -c "import numpy; print(numpy.__version__)"
```

### Dockerfile
Building docker image.
```
docker build -t jakublala/openmm-md:latest .
docker push jakublala/openmm-md:latest
```
Use it in a jobscript with PBS job scheduler.
Note that the current (latest) image supports CUDA 11.8! I am not sure if that has to be the same as the one loaded? Maybe not... :/
```
module load apptainer
apptainer pull docker://jakublala/openmm-md:latest
```

Run interactively (debugging and development):
```
docker run -it jakublala/openmm-md
```


Check some python code in Apptainer
```
apptainer exec openmm-md_v1.0.0.sif python -c "import numpy; print(numpy.__version__)"
```

