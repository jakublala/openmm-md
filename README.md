# openmm-md


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


### Required packages
- `openmm cudatoolkit=11.8`
- `fire`
- `mdtraj`
- `pdbfixer`


- `netCDF4` - speeds up MDanalysis logging of the trajectory


### Getting data from HX1 to Workstation
`scp -r jl24018@login.hx1.hpc.ic.ac.uk:/gpfs/home/jl24018/projects/openmm-md/scripts/241010_FoldingUponBinding/output/output.zip output.zip`