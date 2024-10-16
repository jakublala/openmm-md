from openmm import *
from openmm.app import *
from openmm.unit import *
from openmm import unit

import random
import numpy as np


def opes(filename, mdtime, timestep, device_index='0', temperature=300, restart_checkpoint=None):

    pdf = PDBFile(f'tmp/{filename}/{filename}_solvated.pdb')

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # System Configuration
    nonbondedMethod = PME
    nonbondedCutoff = 1.0*nanometers  #This cutoff and the following for tolerance are standard values
    ewaldErrorTolerance = 10**(-4)
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 10**(-4) #Quite a default value for accuracy

    # Integration Options
    dt = 0.001 * picoseconds * timestep
    temperature = temperature*kelvin
    friction = 1.0/picosecond

    # Simulation Options
    steps = int(mdtime * nanoseconds / dt)
    equilibrationSteps = int(1 * nanosecond / dt)
    platform = Platform.getPlatformByName('CUDA')

    traj_interval = int(100 * picoseconds / dt)

    from mdareporter import MDAReporter
    trajReporter = MDAReporter(
            f'tmp/{filename}/{filename}.dcd', 
            traj_interval, 
            enforcePeriodicBox=False, 
            selection="protein"
            )

    dataReporter = StateDataReporter(
                file=f'tmp/{filename}/{filename}.out',
                reportInterval=traj_interval,
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
                totalSteps=steps,
                separator='\t'
            )
    checkpointReporter = CheckpointReporter(f'tmp/{filename}/{filename}.chk', traj_interval)


    # Prepare the Simulation
    print('Building system...')
    topology = pdf.topology
    positions = pdf.positions
    system = forcefield.createSystem(
        topology, 
        nonbondedMethod=nonbondedMethod, 
        nonbondedCutoff=nonbondedCutoff,
        constraints=constraints, 
        rigidWater=rigidWater, 
        ewaldErrorTolerance=ewaldErrorTolerance
        )

    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    properties = {'DeviceIndex': device_index}
    simulation = Simulation(topology, system, integrator, platform, properties)
    
    if restart_checkpoint: 
        simulation.loadCheckpoint(restart_checkpoint)
        # no equilibration for system from checkpoint
    else:
        simulation.context.setPositions(positions)
    
        #Equilibrate
        print('Equilibrating...')
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(equilibrationSteps)

    # Add PLUMED bias
    with open(f'tmp/{filename}/{filename}_plumed.dat', 'r') as file:
        script = file.read()

    from openmmplumed import PlumedForce
    plumed_force = PlumedForce(script)
    system.addForce(plumed_force)
    
    # Reinitialize the simulation with the updated system
    simulation.system = system
    simulation.context.reinitialize(preserveState=True)

    # Simulate
    print('Simulating...')

    simulation.reporters.append(trajReporter)
    simulation.reporters.append(dataReporter)
    simulation.reporters.append(checkpointReporter)
    simulation.currentStep = 0
    simulation.step(steps)
    

    # TODO: DO ADAPTIVE CONVERGENCE HERE!!!!!
    # i.e. extend the simulation if convergence is not reached


    # move from tmp/ to output/
    import shutil, os
    os.makedirs('output', exist_ok=True)
    shutil.copytree(f'tmp/{filename}', f'output/{filename}', dirs_exist_ok=True)

