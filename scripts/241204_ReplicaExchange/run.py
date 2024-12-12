
from openmm import unit
from openmmtools import testsystems, states, mcmc
from openmmtools.multistate import ParallelTemperingSampler, MultiStateReporter
import tempfile
import os

from openmm.app import PDBFile, ForceField, PME, HBonds

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


def replicate_plumed_file(file_path, n_replicas):
    with open(file_path, 'r') as file:
        plumed_script = file.read()
    for i in range(n_replicas):
        # Replace X with replica number in both hills and colvar filenames
        modified_script = plumed_script.replace('sumo1c_X', f'results/sumo1c_{i}')
        
        with open(f'plumed_{i}.dat', 'w') as file:
            file.write(modified_script)

def main():


    from openmmtools import cache
    from openmm.openmm import Platform
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': str(rank), 'Precision': 'mixed'}
    cache.global_context_cache.set_platform(platform, properties)
    
    # pdb_path = '../../data/241010_FoldingUponBinding/output/SUMO-1C/241128-MetaD/sumo1c_equilibrated.pdb'
    # pdb = PDBFile(pdb_path)

    from openmm.app import PDBxFile
    pdb_path = '../241211_ReplicaImplementation/test/CD28_general_equilibrated.cif'
    pdb = PDBxFile(pdb_path)

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')


    n_replicas = 4  # Number of temperature replicas.
    T_min = 300.0 * unit.kelvin  # Minimum temperature.
    T_max = 600.0 * unit.kelvin  # Maximum temperature.

    timestep = 2.0*unit.femtoseconds
    md_steps = 500 # as suggested in some paper


    PLUMED = True

    from openmmplumed import PlumedForce

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
            replicate_plumed_file('plumed.dat', n_replicas)
            with open(f'plumed_{i}.dat', 'r') as file:
                plumed_script = file.read()
            _system.addForce(PlumedForce(plumed_script))
        systems.append(_system)

    move = mcmc.LangevinDynamicsMove(
        timestep=timestep,
        collision_rate=1/unit.picoseconds,
        n_steps=md_steps # each move is 1 ps
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
    for system in systems:
        sampler_states.append(states.SamplerState(
            positions=pdb.positions,
            box_vectors=system.getDefaultPeriodicBoxVectors()
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
    temperatures = [T_min + (T_max - T_min) * (np.exp(float(i) / float(n_replicas-1)) - 1.0) / (np.e - 1.0) for i in range(n_replicas)]
    logger.info(f"Assigning temperatures: {temperatures}")

    # thermodynamic state holds information on ensemble parameters like temperture and pressure
    thermodynamic_states = [states.ThermodynamicState(system=system, temperature=T) for system, T in zip(systems, temperatures)]

    # with no mpi4py, we are getting a single GPU performance


    simulation.create(
        thermodynamic_states=thermodynamic_states,
        sampler_states=sampler_states, # can be a single state, which gets copied to all replicas
        storage=reporter
        )
    
    simulation.run()

    # TODO: does the swap take into account the PLUMED ENERGY?!?!?!


    # logger.info(f'Accepted steps: {move.n_accepted}')
    # logger.info(f'Proposed steps: {move.n_proposed}')
    # logger.info(f'Fraction accepted: {move.fraction_accepted}')


    # print(analyzer.get_effective_energy_timeseries())




if __name__ == "__main__":
    main()
