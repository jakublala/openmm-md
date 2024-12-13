
from openmm import unit
from openmmtools import testsystems, states, mcmc
from openmmtools.multistate import ParallelTemperingSampler, MultiStateReporter
import tempfile
import os

from openmmplumed import PlumedForce
import numpy as np
import pathlib
from openmmtools.multistate import ReplicaExchangeSampler
import MDAnalysis as mda

from openmmtools import cache
from openmm.openmm import Platform
from openmm.app import PDBxFile

from openmm.app import ForceField, PME, HBonds

import logging
import sys

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


def replicate_plumed_file(output_dir, filename, n_replicas):
    for i in range(n_replicas):
        with open(f'{output_dir}/{filename}_plumed.dat', 'r') as file:
            plumed_script = file.read()

        modified_script = plumed_script.replace(
            f"FILE={output_dir}/{filename}", 
            f"FILE={output_dir}/{filename}_{i}"
            )
            
        with open(f'{output_dir}/{filename}_plumed_{i}.dat', 'w') as file:
            file.write(modified_script)


def run_replica_plumed(
        filename: str,
        mdtime: int,
        timestep: float, 
        temperatures: list[float],
        swap_time: float,
        device: str,
        device_precision: str,
        output_dir: str, 
        logging_frequency: int,
        plumed_config: str,
        chain_mode: str,
        n_replicas: int,
        ):
    
    from src.plumed.io import create_plumed_input
    create_plumed_input(
        filename=filename,
        output_dir=output_dir,
        config=plumed_config,
        mode=chain_mode,
    )
    replicate_plumed_file(output_dir, filename, n_replicas)
    
    if device != 'cuda':
        raise NotImplementedError("Replica exchange only supports CUDA")
    
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': str(rank), 'Precision': device_precision}
    cache.global_context_cache.set_platform(platform, properties)
    
    pdb = PDBxFile(f"{output_dir}/{filename}_equilibrated.cif")

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    n_replicas = len(temperatures)  # Number of temperature replicas.

    timestep = timestep * unit.femtoseconds
    swap_time = swap_time * unit.picoseconds
    swap_steps = int(swap_time / timestep)

    md_time = mdtime * unit.nanosecond
    md_steps = int(md_time / timestep)
    N_ITERATIONS = int(md_steps / swap_steps)

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
            with open(f'plumed_{i}.dat', 'r') as file:
                plumed_script = file.read()
            _system.addForce(PlumedForce(plumed_script))
        systems.append(_system)

    move = mcmc.LangevinDynamicsMove(
        timestep=timestep,
        collision_rate=1/unit.picoseconds,
        n_steps=md_steps
    )

    storage_path = pathlib.Path(f"{output_dir}/replica_exchange.nc")
    if rank == 0:
        # if exists, delete
        if storage_path.exists():
            storage_path.unlink()
        checkpoint_path = pathlib.Path(f"{output_dir}/replica_exchange_checkpoint.nc")
        if checkpoint_path.exists():
            checkpoint_path.unlink()
    
    # Instantiate the reporter that only collects protein atoms
    reporter = MultiStateReporter(
        storage=storage_path, 
        checkpoint_interval=10,
        analysis_particle_indices=mda.Universe(pdb).select_atoms("protein").ids
        )

    # Create sampler state with both positions and box vectors
    # this contains the State, meaning the positions, velocities, and box vectors, etc.
    sampler_states = []
    for system in systems:
        sampler_states.append(states.SamplerState(
            positions=pdb.positions,
            box_vectors=system.getDefaultPeriodicBoxVectors()
            # hopefully velocities if None are set then by the temperature later
            # TODO: check if true
        ))

    simulation = ReplicaExchangeSampler(
        replica_mixing_scheme='swap-all',
        # or swap-all, which is more expensive
        # and does n_replicas**3 swap attempts per iteration
        mcmc_moves=move, 
        number_of_iterations=N_ITERATIONS # this does 50 iterations of the move, hence 50 * 500 steps
        )
    
    # thermodynamic state holds information on ensemble parameters like temperture and pressure
    thermodynamic_states = [states.ThermodynamicState(system=system, temperature=T) for system, T in zip(systems, temperatures)]

    # with no mpi4py, we are getting a single GPU performance


    simulation.create(
        thermodynamic_states=thermodynamic_states,
        sampler_states=sampler_states, # can be a single state, which gets copied to all replicas
        storage=reporter
        )
    
    simulation.run()
