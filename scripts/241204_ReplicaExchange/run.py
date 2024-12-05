
from openmm import unit
from openmmtools import testsystems, states, mcmc
from openmmtools.multistate import ParallelTemperingSampler, MultiStateReporter
import tempfile
import os

from openmm.app import PDBFile, ForceField, PME, HBonds

import logging
import sys
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
# Set up basic logging to console
logging.basicConfig(
    level=logging.INFO,  # Change to DEBUG for more detailed output
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)

def main():


    from openmmtools import cache
    from openmm.openmm import Platform
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': '0, 1', 'Precision': 'mixed'}
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
    # if exists, delete
    if storage_path.exists():
        storage_path.unlink()
    checkpoint_path = pathlib.Path('replica_exchange_checkpoint.nc')
    if checkpoint_path.exists():
        checkpoint_path.unlink()
    
    reporter = MultiStateReporter(
        storage=storage_path, 
        checkpoint_interval=10
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

    from openmmtools.multistate import ParallelTemperingAnalyzer

    # analyzer = ParallelTemperingAnalyzer(reporter)
    # print(analyzer.read_energies())

    # open .nc file?!?!
    import netCDF4
    nc = netCDF4.Dataset(storage_path, 'r')
    print(nc.variables.keys())

    # loading traj and write to dcd format
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader
    coord_traj = np.zeros((n_replicas, 20, pdb.topology.getNumAtoms(), 3)) * unit.nanometer
    for time in range(0, N_ITERATIONS, 10):
        state_samplers = reporter.read_sampler_states(time)
        for i, state_sampler in enumerate(state_samplers):
            coord_traj[i, int(time/10)] += state_sampler.positions

    # only taking the first replica trajectory
    u = mda.Universe(pdb_path, coord_traj[0].value_in_unit(unit.angstrom), format=MemoryReader)
    with mda.Writer("test.dcd", u.select_atoms("protein").n_atoms) as W:
        for ts in u.trajectory:
            W.write(u.select_atoms("protein"))

    # print(analyzer.get_effective_energy_timeseries())




if __name__ == "__main__":
    main()
