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


rank, n_procs = _setup_mpi()
logger = _setup_logging(rank)


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

    _ = _get_platform(device, rank)
    # TIMESTEP = timestep * unit.femtoseconds
    TIMESTEP = 2.0 * unit.femtoseconds
    SWAP_TIME = swap_time * unit.picoseconds
    SWAP_STEPS = int(SWAP_TIME / TIMESTEP)
    N_STEPS = int(mdtime / timestep)
    N_ITERATIONS = int(N_STEPS / swap_time)

    cif = PDBxFile(get_file_by_extension(output_dir, '_equilibrated.cif'))

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    N_REPLICAS = len(temperatures)
    logger.info(f"Running {N_REPLICAS} replicas with temperatures {temperatures}")

    # assert temperatures elements have unit of unit.kelvin
    assert all(isinstance(temp, unit.Quantity) and temp.unit == unit.kelvin for temp in temperatures), "Temperatures must have units of kelvin"
    logger.info(f"Assigning temperatures: {temperatures}")

    # System Configuration
    nonbondedMethod = PME
    nonbondedCutoff = 1.0*unit.nanometers  # This cutoff and the following for tolerance are standard values
    ewaldErrorTolerance = 10**(-4)
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 10**(-4)

    # Create systems and add independent plumed forces to each
    systems = []
    for i in range(N_REPLICAS):
        _system = forcefield.createSystem(
            topology=cif.topology, 
            nonbondedMethod=nonbondedMethod, 
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints, 
            rigidWater=rigidWater, 
            ewaldErrorTolerance=ewaldErrorTolerance,
        )
        if rank == 0:
            plumed_filepath = get_file_by_extension(output_dir, f'{filename}_plumed.dat')
            _replicate_plumed_file(filename, plumed_filepath, output_dir, i)
            
        with open(f'{output_dir}/{filename}_plumed_{i}.dat', 'r') as file:
            plumed_script = file.read()
        # _system.addForce(PlumedForce(plumed_script))
        systems.append(_system)
    
    
    # if rank == 0:
    #     os.remove(plumed_filepath)

    # we don't use this, but let's keep it for reference
    # move = mcmc.GHMCMove(
    #     timestep=timestep,  # Each integration step is 2 fs (Langevin dynamics)
    #     n_steps=md_steps,
    #     collision_rate=1/unit.picoseconds, # equivalent to friction coefficient in LangevinMiddleIntegrator
    # )    

    # Define Monte Carlo move
    move = mcmc.LangevinDynamicsMove(
        timestep=timestep,
        collision_rate=1/unit.picoseconds,
        n_steps=SWAP_STEPS,
        constraint_tolerance=constraintTolerance
    )

    thermodynamic_states = [
        states.ThermodynamicState(system=system, temperature=T) 
        for system, T in zip(systems, temperatures)
    ]
    # Create sampler state with both positions and box vectors
    # this contains the State, meaning the positions, velocities, and box vectors, etc.
    sampler_states = []
    for i, system in enumerate(systems):
        sampler_state = states.SamplerState(
            positions=cif.positions,
            box_vectors=system.getDefaultPeriodicBoxVectors(),
        )

        print(f"Kinetic energy: {sampler_state.kinetic_energy}")
        print(f"Potential energy: {sampler_state.potential_energy}")
        sampler_states.append(sampler_state)


    simulation = ReplicaExchangeSampler(
        replica_mixing_scheme='swap-all',
        # or swap-all, which is more expensive
        # and does n_replicas**3 swap attempts per iteration
        mcmc_moves=move, 
        number_of_iterations=N_ITERATIONS # this does 50 iterations of the move, hence 50 * 500 steps
        )


    storage_path = pathlib.Path(f'{output_dir}/{filename}_replica_exchange.nc')
    _nc_cleanup(rank, storage_path)

    # Use MPI barrier before creating reporter to ensure cleanup is done
    MPI.COMM_WORLD.Barrier()

    reporter = MultiStateReporter(
        storage=storage_path, 
        checkpoint_interval=10,
        analysis_particle_indices=mda.Universe(cif).select_atoms('protein').ids
        )
    
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










