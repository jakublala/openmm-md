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