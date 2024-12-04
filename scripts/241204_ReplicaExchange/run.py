from openmm import unit
from openmmtools import testsystems, states, mcmc
from openmmtools.multistate import ParallelTemperingSampler, MultiStateReporter
import tempfile
import os

def main():
    pdb = PDBFile('/home/jakub/phd/openmm-md/data/241010_FoldingUponBinding/output/SUMO-1C/241128-MetaD/sumo1c_equilibrated.pdb')

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

    move = mcmc.GHMCMove(timestep=2.0*unit.femtoseconds, n_steps=50)


    storage_path = 'replica_exchange.nc'
    reporter = MultiStateReporter(storage_path, checkpoint_interval=10)


    # Create sampler state with both positions and box vectors
    sampler_state = states.SamplerState(
        positions=pdb.positions,
        box_vectors=system.getDefaultPeriodicBoxVectors()
    )

    simulation = ParallelTemperingSampler(
        mcmc_moves=move, 
        number_of_iterations=2
        )

    simulation.create(
        thermodynamic_state=reference_state,
        sampler_states=[sampler_state],
        storage=reporter, 
        min_temperature=T_min,
        max_temperature=T_max, 
        n_temperatures=n_replicas)


    simulation.run()


if __name__ == "__main__":
    main()
