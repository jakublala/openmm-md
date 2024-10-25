# openmm-md

### Install
Create env.
```
conda create -n openmm
conda install python pdbfixer
conda install -c conda-forge openmm
pip install -e .
```

Run `pip install -e .`
Then also run `conda install cudatoolkit=11.8` based on your CUDA version to speed up OpenMM.


```
git clone https://github.com/plumed/plumed2
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