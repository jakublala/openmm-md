
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
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    n_procs = comm.Get_size()
except ImportError:
    rank = 0
    n_procs = 1

# Configure logging only for rank 0
if rank == 0:
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        stream=sys.stdout
    )

def main():


    from openmmtools import cache
    from openmm.openmm import Platform
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': str(rank), 'Precision': 'mixed'}
    cache.global_context_cache.set_platform(platform, properties)
    
    pdb_path = '../../data/241010_FoldingUponBinding/output/SUMO-1C/241128-MetaD/sumo1c_equilibrated.pdb'

    pdb = PDBFile(pdb_path)

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    system = forcefield.createSystem(
        topology=pdb.topology, 
        nonbondedMethod=PME,
        nonbondedCutoff=1*unit.nanometer,
        constraints=HBonds
        )
    

    n_replicas = 4  # Number of temperature replicas.
    T_min = 300.0 * unit.kelvin  # Minimum temperature.
    T_max = 600.0 * unit.kelvin  # Maximum temperature.


    # thermodynamic state holds information on ensemble parameters like temperture and pressure
    reference_state = states.ThermodynamicState(system=system, temperature=T_min)

    timestep = 2.0*unit.femtoseconds
    md_steps = 500 # as suggested in some paper

    move = mcmc.GHMCMove(
        timestep=timestep,  # Each integration step is 2 fs (Langevin dynamics)
        n_steps=md_steps,
        collision_rate=1/unit.picoseconds, # equivalent to friction coefficient in LangevinMiddleIntegrator
    )

    import pathlib
    storage_path = pathlib.Path('replica_exchange.nc')
    if rank == 0:
        # if exists, delete
        if storage_path.exists():
            storage_path.unlink()
        checkpoint_path = pathlib.Path('replica_exchange_checkpoint.nc')
        if checkpoint_path.exists():
            checkpoint_path.unlink()
    

    # get atoms of all protein
    import MDAnalysis as mda
    u = mda.Universe(pdb_path)
    protein_atom_ids = u.select_atoms("protein").ids

    reporter = MultiStateReporter(
        storage=storage_path, 
        checkpoint_interval=10,
        analysis_particle_indices=protein_atom_ids
        )

    # Create sampler state with both positions and box vectors
    sampler_state = states.SamplerState(
        positions=pdb.positions,
        box_vectors=system.getDefaultPeriodicBoxVectors()
        # hopefully velocities if None are set then by the temperature later
        # TODO: check if true
    )
    
    N_ITERATIONS = 10
    simulation = ParallelTemperingSampler(
        mcmc_moves=move, 
        number_of_iterations=N_ITERATIONS # this does 50 iterations of the move, hence 50 * 500 steps
        )
    
    # with no mpi4py, we are getting a single GPU performance

    simulation.create(
        thermodynamic_state=reference_state,
        sampler_states=sampler_state, # can be a single state, which gets copied to all replicas
        storage=reporter, 
        min_temperature=T_min,
        max_temperature=T_max, 
        n_temperatures=n_replicas)
    
    simulation.run()


    # print(analyzer.get_effective_energy_timeseries())




if __name__ == "__main__":
    main()
