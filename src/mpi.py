from mpi4py import MPI
import os
from typing import Optional, Union, List, Tuple
from src.analysis.utils import get_file_by_extension

class MPIContext:
    def __init__(self):
        self.comm = None
        self.rank = 0
        self.n_procs = 1
        
        try:
            # Check if MPI is already initialized
            if not MPI.Is_initialized():
                # Initialize MPI with thread support
                required = MPI.THREAD_MULTIPLE
                provided = MPI.Init_thread(required)
                if provided < required:
                    print(f"Warning: MPI thread support level {provided} is less than required {required}")
            
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.n_procs = self.comm.Get_size()
        except ImportError:
            self.rank = 0
            self.n_procs = 1

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if MPI.Is_initialized():
            MPI.Finalize()