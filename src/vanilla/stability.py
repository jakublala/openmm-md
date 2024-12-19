from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin, nanosecond
import fire
import time
from mdareporter import MDAReporter
from openmm.app import PDBxFile
from src.utils import get_checkpoint_interval, get_platform_and_properties


# Define a function to check if a residue is a water molecule or an ion
def is_water_or_ion(residue):
    if residue.name == 'HOH' or residue.name in ['NA', 'CL']:
        return True
    return False

import logging
import os
logger = logging.getLogger(__name__) 


def stability(
        filename=None,
        device_index='0',
        mdtime=None, # give in ns
        timestep=2, # in femtoseconds
        temperature=300,
        restart=False,
        logging_frequency=100,
        device="cuda",
        output_dir=None,
        equilibrate_only=False,
        ):
    
    dt = 0.001 * timestep * picosecond

    if output_dir is None:
        raise ValueError('Output directory is required')
    if filename is None:
        raise ValueError('Filename is required')
    
    nsteps = int(mdtime * nanosecond / dt)
    
    if restart:
        raise ValueError('Restart not implemented yet')
    else:

        if os.path.exists(f'{output_dir}/{filename}_equilibrated.cif'):
            logger.info(f'Equilibrated state found at {output_dir}/{filename}_equilibrated.cif')
            logger.info('Skipping equilibration...')
            pdb = PDBxFile(f'{output_dir}/{filename}_equilibrated.cif')
            equilibrated = True
        else:
            pdb = PDBxFile(f'{output_dir}/{filename}_solvated.cif')
            equilibrated = False


    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # System Configuration
    nonbondedMethod = PME
    nonbondedCutoff = 1.0*nanometer  # This cutoff and the following for tolerance are standard values
    ewaldErrorTolerance = 10**(-4)
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 10**(-4) # quite a default value for accuracy

    # Integration Options
    dt = 0.001 * picoseconds * timestep
    temperature = temperature * kelvin
    friction = 1.0/picoseconds

    system = forcefield.createSystem(
        topology=pdb.topology, 
        nonbondedMethod=nonbondedMethod,
        nonbondedCutoff=nonbondedCutoff,
        constraints=constraints,
        rigidWater=rigidWater,
        ewaldErrorTolerance=ewaldErrorTolerance
        )
    print('Periodic boundary conditions:', system.usesPeriodicBoundaryConditions())

    # TODO: run it at body temperature
    # have a bit larger friction coefficient initially in an equilibriation phase
    # then decrease it to a smaller value for production run
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)

    platform, properties = get_platform_and_properties(device, device_index)

    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    

    positions = pdb.positions
    equilibrationSteps = int(1 * nanosecond / dt)
    if restart:
        raise ValueError('Restart not implemented yet')
    else:
        simulation.context.setPositions(pdb.positions)

    if equilibrated:
        if equilibrate_only:
            raise ValueError("Equilibrate_only is True, but equilibrated state found. You don't need to equilibrate again.")
        simulation.context.setPositions(positions)
    else:
        simulation.context.setPositions(positions)

        # Equilibrate
        logger.info('Equilibrating...')
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(equilibrationSteps)

        from src.utils import save_equilibrated_state
        save_equilibrated_state(
            simulation=simulation,
            output_dir=output_dir,
            filename=filename
            )
        
        if equilibrate_only:
            return
    
    simulation.reporters.append(
        MDAReporter(
            f'{output_dir}/{filename}.dcd', 
            logging_frequency, 
            enforcePeriodicBox=False, 
            selection="protein"
            )
        )
    
    # log the energy and temperature every 1000 steps
    simulation.reporters.append(
        StateDataReporter(
            file=f'{output_dir}/{filename}.out',
            reportInterval=logging_frequency,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            progress=True,
            remainingTime=True,
            totalSteps=nsteps
        )
    )

    simulation.reporters.append(
        CheckpointReporter(
            file=f'{output_dir}/{filename}.chk', 
            reportInterval=get_checkpoint_interval(timestep)
            )
        )

    simulation.step(nsteps)

if __name__ == '__main__':
    fire.Fire(stability)