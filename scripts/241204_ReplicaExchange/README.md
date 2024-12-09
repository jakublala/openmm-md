### REST



To run with two GPUs
```
export OMP_NUM_THREADS=32
mpirun -np 2 --hostfile hostfile python run.py
```

Haven't tested whether increasing the number of CPU threads does any better