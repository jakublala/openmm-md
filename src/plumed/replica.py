
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
    # Add these lines:
    logging.getLogger('openmmtools').setLevel(logging.DEBUG)
    logging.getLogger('openmmtools.multistate').setLevel(logging.DEBUG)
else:
    class NoOpLogger:
        def __getattr__(self, name):
            def method(*args, **kwargs):
                pass
            return method

    logger = NoOpLogger()


def replicate_plumed_file(output_dir, filename, temperatures):
    for i, T in enumerate(temperatures):
        with open(f'{output_dir}/{filename}_plumed.dat', 'r') as file:
            plumed_script = file.read()

        modified_script = plumed_script.replace(
            f"FILE={output_dir}/{filename}", 
            f"FILE={output_dir}/{filename}_{i}"
            )
        
        modified_script = modified_script.replace(
            f"TEMP=300", 
            f"TEMP={T.value_in_unit(unit.kelvin)}"
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
        replica_type: str,
        ):
    
    from src.plumed.io import create_plumed_input
    create_plumed_input(
        filename=filename,
        output_dir=output_dir,
        config=plumed_config,
        mode=chain_mode,
    )
    replicate_plumed_file(output_dir, filename, temperatures)
    
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

    import openmm as mm

    class NoPeriodicCVLangevinDynamicsMove(mcmc.LangevinDynamicsMove):
        def _after_integration(self,
            context: mm.Context,
            thermodynamic_state: states.ThermodynamicState
            ):
            """Wrap coordinates into primary unit cell after integration."""
            system: mm.System = thermodynamic_state.get_system()
            if system.usesPeriodicBoundaryConditions():
                state = context.getState(getPositions=True, enforcePeriodicBox=True)
                positions = state.getPositions()
                box_vectors = state.getPeriodicBoxVectors()
                
                # we need to move the centre of each of these chains to be in the original, same unit cell
                chain_A_positions = positions[1:504]
                chain_B_positions = positions[505:2714]

                # turn into a numpy array
                chain_A_positions = np.array(chain_A_positions)
                chain_B_positions = np.array(chain_B_positions)

                chain_A_centre = np.mean(chain_A_positions, axis=0)
                chain_B_centre = np.mean(chain_B_positions, axis=0)


                with open(f"{output_dir}/chain_centres.txt", "a") as file:
                    file.write(str(chain_A_centre[0].value_in_unit(unit.nanometer)))
                    file.write("\t")
                    file.write(str(chain_A_centre[1].value_in_unit(unit.nanometer)))
                    file.write("\t")
                    file.write(str(chain_A_centre[2].value_in_unit(unit.nanometer)))
                    file.write("\t")
                    file.write("\t")
                    
                    file.write(str(chain_B_centre[0].value_in_unit(unit.nanometer)))
                    file.write("\t")
                    file.write(str(chain_B_centre[1].value_in_unit(unit.nanometer)))
                    file.write("\t")
                    file.write(str(chain_B_centre[2].value_in_unit(unit.nanometer)))
                    file.write("\n")

            # Call parent method in case future versions add functionality
            super()._after_integration(context, thermodynamic_state)

    if plumed_config['cv1.pbc'] and plumed_config['cv2.pbc']:
        move = mcmc.LangevinDynamicsMove(
            timestep=timestep,
            collision_rate=1/unit.picoseconds,
            n_steps=swap_steps,
            reassign_velocities=True
        )
    else:
        move = NoPeriodicCVLangevinDynamicsMove(
            timestep=timestep,
            collision_rate=1/unit.picoseconds,
            n_steps=swap_steps,
            reassign_velocities=True
        )

    # TODO: adjust this move to move always to base image?

    storage_path = pathlib.Path(f"{output_dir}/replica_exchange.nc")
    checkpoint_path = pathlib.Path(f"{output_dir}/replica_exchange_checkpoint.nc")
    if rank == 0:
        # if exists, delete
        if storage_path.exists():
            storage_path.unlink()
        if checkpoint_path.exists():
            checkpoint_path.unlink()
    
    # Instantiate the reporter that only collects protein atoms
    reporter = MultiStateReporter(
        storage=storage_path,
        checkpoint_storage='replica_exchange_checkpoint.nc',
        checkpoint_interval=100,
        analysis_particle_indices=mda.Universe(pdb).select_atoms("protein").ids
        )

    # Create sampler state with both positions and box vectors
    # this contains the State, meaning the positions, velocities, and box vectors, etc.
    # sampler_states = []
    # import copy
    # for system in systems:
    #     positions = copy.deepcopy(pdb.positions)
    #     sampler_states.append(states.SamplerState(
    #         positions=positions,
    #         box_vectors=system.getDefaultPeriodicBoxVectors()
    #         # hopefully velocities if None are set then by the temperature later
    #         # TODO: check if true
    #     ))

    simulation = ReplicaExchangeSampler(
        replica_mixing_scheme='swap-all',
        # or swap-all, which is more expensive
        # and does n_replicas**3 swap attempts per iteration
        mcmc_moves=move, 
        number_of_iterations=N_ITERATIONS # this does 50 iterations of the move, hence 50 * 500 steps
        )
    

    # get thermodynamic states and sampler states
    if replica_type == 'rest2':
        from src.replica.rest2 import get_thermodynamic_and_sampler_states
        system = forcefield.createSystem(
            topology=pdb.topology, 
            nonbondedMethod=PME,
            nonbondedCutoff=1*unit.nanometer,
            constraints=HBonds
            )
        
        thermodynamic_states, sampler_states = get_thermodynamic_and_sampler_states(
            system=system, 
            temperatures=temperatures,
            equilibriated_file=f'{output_dir}/{filename}_equilibrated.cif',
            plumed_input_file=f'{output_dir}/{filename}_plumed.dat'
            )
    else:
        # system holds information on the force field parameters
        systems = []
        for i in range(n_replicas):
            _system = forcefield.createSystem(
                topology=pdb.topology, 
                nonbondedMethod=PME,
                nonbondedCutoff=1*unit.nanometer,
                constraints=HBonds
                )
        
            with open(f'{output_dir}/{filename}_plumed_{i}.dat', 'r') as file:
                plumed_script = file.read()
            _system.addForce(PlumedForce(plumed_script))
            systems.append(_system)
        # thermodynamic state holds information on ensemble parameters like temperture and pressure
        thermodynamic_states = [
            states.ThermodynamicState(system=system, temperature=T) 
            for system, T in zip(systems, temperatures)
        ]
        sampler_states = [states.SamplerState(
            positions=pdb.positions,
            box_vectors=system.getDefaultPeriodicBoxVectors()
        ) for system in systems]

    # with no mpi4py, we are getting a single GPU performance
    logger.info(f"Running {n_replicas} replicas with temperatures {temperatures}")
    logger.info(f"Running with timestep {timestep} and mdtime {mdtime * unit.nanoseconds}, which is {N_ITERATIONS} iterations of replica swaps.")
    logger.info(f"Replica swaps every {swap_time}, which is {swap_steps} steps")



    simulation.create(
        thermodynamic_states=thermodynamic_states,
        sampler_states=sampler_states, # can be a single state, which gets copied to all replicas
        storage=reporter
        )
            
    # # Equilibriate the replicas at the new temperatures
    # equilibriation_time = 10 * unit.nanoseconds
    # equilibriation_steps = int(equilibriation_time / swap_time)
    # simulation.equilibrate(n_iterations=equilibriation_steps, mcmc_moves=move)
    
    simulation.run()
