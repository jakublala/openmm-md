from openmm import unit
from openmmtools import testsystems, states, mcmc
from openmmtools.multistate import ReplicaExchangeSampler
from openmmtools.multistate import ParallelTemperingSampler, MultiStateReporter
from openmmplumed import PlumedForce

from typing import List
import pathlib
import tempfile
import os
from src.mpi import MPIContext

from src.analysis.utils import get_file_by_extension
from src.plumed.io import create_plumed_input

from openmm.app import PDBxFile, ForceField, PME, HBonds

from openmmtools import cache
from openmm.openmm import Platform

import logging
import sys
import numpy as np
from mpi4py import MPI
import MDAnalysis as mda

# Set up MPI
try:
    from mpi4py import MPI
    # Check if MPI is already initialized
    if not MPI.Is_initialized():
        # Initialize MPI with thread support
        required = MPI.THREAD_MULTIPLE
        provided = MPI.Init_thread(required)
        if provided < required:
            print(f"Warning: MPI thread support level {provided} is less than required {required}")
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    n_procs = comm.Get_size()
except ImportError:
    rank = 0
    n_procs = 1

# Configure logging
if rank == 0:
    logger = logging.getLogger()
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        stream=sys.stdout
    )
else:
    class NoOpLogger:
        def __getattr__(self, name):
            def method(*args, **kwargs):
                pass
            return method

    logger = NoOpLogger()

# TODO: for now, keep it here, but then move elsewhere to MPI utils
def _setup_mpi():
    # Set up MPI
    try:
        # Check if MPI is already initialized
        if not MPI.Is_initialized():
            # Initialize MPI with thread support
            required = MPI.THREAD_MULTIPLE
            provided = MPI.Init_thread(required)
            if provided < required:
                print(f"Warning: MPI thread support level {provided} is less than required {required}")
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        n_procs = comm.Get_size()
    except ImportError:
        rank = 0
        n_procs = 1
    return rank, n_procs

def _setup_logging(rank):
    # Configure logging
    logger = logging.getLogger(__name__)  # Get logger for this module specifically
    if rank == 0:
        # Clear any existing handlers to avoid duplicate logging
        logger.handlers.clear()
        
        # Create console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        console_handler.setFormatter(formatter)
        
        # Add handler to logger
        logger.addHandler(console_handler)
        logger.setLevel(logging.DEBUG)
        
        # Prevent propagation to root logger
        logger.propagate = False
    else:
        class NoOpLogger:
            def __getattr__(self, name):
                def method(*args, **kwargs):
                    pass
                return method
        logger = NoOpLogger()
    return logger


def _get_platform(device, rank):
    assert device == 'cuda', "Only CUDA is supported for Replica Exchange for now"
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': str(rank), 'Precision': 'mixed'}
    cache.global_context_cache.set_platform(platform, properties)
    return platform

def _replicate_plumed_file(filename, file_path, output_dir, replica_id):
    with open(file_path, 'r') as file:
        plumed_script = file.read()

    modified_script = plumed_script.replace(
        f"FILE={output_dir}/{filename}", 
        f"FILE={output_dir}/{filename}_{replica_id}"
        )
        
    with open(f'{output_dir}/{filename}_plumed_{replica_id}.dat', 'w') as file:
            file.write(modified_script)

def _nc_cleanup(rank, storage_path):
    if rank == 0:
        # if exists, delete
        if storage_path.exists():
            logger.info(f"Deleting {storage_path} for next replica exchange simulation")
            storage_path.unlink()
        checkpoint_path = pathlib.Path(str(storage_path).replace('.nc', '') + '_checkpoint.nc')
        if checkpoint_path.exists():
            checkpoint_path.unlink()
    
def run_replica_plumed(
        filename: str,
        mdtime: int,
        timestep: float,
        swap_time: int,
        temperatures: List[float],
        device: str,
        device_precision: str,
        output_dir: str,
        logging_frequency: int,
        plumed_config: dict,
        chain_mode: str,
        ):
    

    create_plumed_input(
        filename=filename, 
        output_dir=output_dir,
        config=plumed_config,
        mode=chain_mode,
        )

    logger.info(f"Running replica {rank} of {n_procs} replicas")

    # _ = _get_platform(device, rank)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': str(rank), 'Precision': 'mixed'}
    cache.global_context_cache.set_platform(platform, properties)
    
    # TIMESTEP = timestep * unit.femtoseconds
    TIMESTEP = 0.01 * unit.femtoseconds
    SWAP_TIME = swap_time * unit.picoseconds
    SWAP_STEPS = int(SWAP_TIME / TIMESTEP)
    N_STEPS = int(mdtime / timestep)
    N_ITERATIONS = int(N_STEPS / swap_time)

    cif = PDBxFile(get_file_by_extension(output_dir, '_equilibrated.cif'))

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # N_REPLICAS = len(temperatures)
    # logger.info(f"Running {N_REPLICAS} replicas with temperatures {temperatures}")

    # assert temperatures elements have unit of unit.kelvin
    # assert all(isinstance(temp, unit.Quantity) and temp.unit == unit.kelvin for temp in temperatures), "Temperatures must have units of kelvin"
    logger.info(f"Assigning temperatures: {temperatures}")


    from openmmplumed import PlumedForce

    pdb = PDBxFile(get_file_by_extension(output_dir, '_equilibrated.cif'))
    n_replicas = 4

    PLUMED = False

    # system holds information on the force field parameters
    systems = []
    for i in range(n_replicas):
        _system = forcefield.createSystem(
            topology=pdb.topology, 
            nonbondedMethod=PME,
            nonbondedCutoff=1*unit.nanometer,
            constraints=HBonds
            )
        if PLUMED:
            plumed_filepath = f'{output_dir}/{filename}_plumed.dat'
            _replicate_plumed_file(f'{output_dir}/plumed.dat', output_dir, i)
            with open(f'{output_dir}/plumed_{i}.dat', 'r') as file:
                plumed_script = file.read()
            _system.addForce(PlumedForce(plumed_script))
        systems.append(_system)

    move = mcmc.LangevinDynamicsMove(
        timestep=timestep,
        collision_rate=1/unit.picoseconds,
        n_steps=500, # each move is 1 ps
        reassign_velocities=True
    )

    import pathlib
    storage_path = pathlib.Path('results/replica_exchange.nc')
    if rank == 0:
        # if exists, delete
        if storage_path.exists():
            storage_path.unlink()
        checkpoint_path = pathlib.Path('results/replica_exchange_checkpoint.nc')
        if checkpoint_path.exists():
            checkpoint_path.unlink()
    

    # get atoms of all protein
    import MDAnalysis as mda
    u = mda.Universe(pdb)
    protein_atom_ids = u.select_atoms("protein").ids

    reporter = MultiStateReporter(
        storage=storage_path, 
        checkpoint_interval=10,
        analysis_particle_indices=protein_atom_ids
        )

    # Create sampler state with both positions and box vectors
    # this contains the State, meaning the positions, velocities, and box vectors, etc.
    sampler_states = []
    import numpy as np
    for system in systems:
        sampler_states.append(states.SamplerState(
            positions=pdb.positions,
            box_vectors=system.getDefaultPeriodicBoxVectors(),
            velocities=np.zeros((system.getNumParticles(), 3)) * unit.nanometers / unit.picoseconds
            # hopefully velocities if None are set then by the temperature later
            # TODO: check if true
        ))

    # run for 1 ns = 1000 ps
    N_ITERATIONS = 1000
    from openmmtools.multistate import ReplicaExchangeSampler
    simulation = ReplicaExchangeSampler(
        replica_mixing_scheme='swap-all',
        # or swap-all, which is more expensive
        # and does n_replicas**3 swap attempts per iteration
        mcmc_moves=move, 
        number_of_iterations=N_ITERATIONS # this does 50 iterations of the move, hence 50 * 500 steps
        )
    
    import numpy as np
    T_min = 300.0 * unit.kelvin  # Minimum temperature.
    T_max = 600.0 * unit.kelvin  # Maximum temperature.

    temperatures = [T_min + (T_max - T_min) * (np.exp(float(i) / float(n_replicas-1)) - 1.0) / (np.e - 1.0) for i in range(n_replicas)]
    logger.info(f"Assigning temperatures: {temperatures}")

    # thermodynamic state holds information on ensemble parameters like temperture and pressure
    thermodynamic_states = [
        states.ThermodynamicState(system=system, temperature=T) 
        for system, T in zip(systems, temperatures)
        ]

    # with no mpi4py, we are getting a single GPU performance


    simulation.create(
        thermodynamic_states=thermodynamic_states,
        sampler_states=sampler_states, # can be a single state, which gets copied to all replicas
        storage=reporter
        )
    
    simulation.run()

def run_replica_task(replica_id, *args, **kwargs):
    """
    Function to handle a single replica's task.
    """
    print(f"Running replica {replica_id} at temperature {kwargs['temperature']}")
    # Call your original replica-running function here with the appropriate parameters
    # Example:
    run_replica_plumed(*args, **kwargs)
    print(f"Replica {replica_id} completed.")
    return f"Replica {replica_id} finished at temperature {kwargs['temperature']}"




if __name__ == '__main__':
    config = {
        'type': 'metad',
        'temperature': 300,
        'stride': 1,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': 1,
        'metad.pace': 500,
        'cv1.type': 'cmap',
        'cv1.sigma': 0.15,
        'cv1.grid_min': None,
        'cv1.grid_max': None,
        'cv1.grid_bin': None,
        # 'cv1.grid_min': 0,
        # 'cv1.grid_max': 45,
        # 'cv1.grid_bin': 200,
        'cv1.pbc': False,
        'cv2.type': 'd',
        'cv2.sigma': 0.27,
        'cv2.grid_min': None,
        'cv2.grid_max': None,
        'cv2.grid_bin': None,
        # 'cv2.grid_min': 0,
        # 'cv2.grid_max': 12,
        # 'cv2.grid_bin': 200,
        'cv2.pbc': False,
        'metad.height': 1.25, # 1/2 * kBT
        'metad.biasfactor': 48,
        'upper_wall.at': 6, # keep this at UW=5, we are primarily looking at BIASFACTOR now
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': None,
        'spot2_residues': None,
        'idr_residues': None,
        'restart': False,
        'trajectory_logging': True
    }

    run_replica_plumed(
        filename='CD28_general',
        mdtime=1000,
        timestep=2,
        swap_time=1,
        temperatures=None,
        device='cuda',
        device_precision='mixed',
        output_dir='/app/scripts/241211_ReplicaImplementation/test',
        logging_frequency=1,
        plumed_config=config,
        chain_mode='two-chain',
    )

