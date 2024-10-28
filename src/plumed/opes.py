from openmm import *
from openmm.app import *
from openmm.unit import *


import logging
logger = logging.getLogger(__name__)


def opes(
        filename, 
        mdtime, 
        timestep, 
        device_index='0', 
        temperature=300, 
        restart_checkpoint=None, 
        device='cuda', 
        output_dir=None, 
        chk_interval=None
        ):

    if output_dir is None:
        raise ValueError('Output directory is required')

    pdf = PDBFile(f'{output_dir}/{filename}_solvated.pdb')

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

    from openmm import Platform

    # List all available platforms
    logger.info("Available Platforms:")
    for i in range(Platform.getNumPlatforms()):
        platform = Platform.getPlatform(i)
        logger.info(f"{i+1}. {platform.getName()}")
    

    if device == 'cuda':
        platform = Platform.getPlatformByName('CUDA')
    elif device == 'cpu':
        platform = Platform.getPlatformByName('CPU')
    else:
        raise ValueError(f'Invalid device: {device}')

    traj_interval = int(100 * picoseconds / dt)

    from mdareporter import MDAReporter
    trajReporter = MDAReporter(
            f'{output_dir}/{filename}.dcd', 
            traj_interval, 
            enforcePeriodicBox=False, 
            selection="protein"
            )

    dataReporter = StateDataReporter(
                file=f'{output_dir}/{filename}.out',
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
    if chk_interval is None:
        logger.warning("Checkpoint interval not specified, using default of every 1 ns")
        chk_interval = int(1 * nanoseconds / dt)
    checkpointReporter = CheckpointReporter(f'{output_dir}/{filename}.chk', chk_interval)


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
    
    if device_index == 'nan':
        simulation = Simulation(topology, system, integrator, platform)
    else:
        logger.info(f"Using device index {device_index} of type {type(device_index)}")
        device_index = "0,1" # HACK: as these seem to be already taking into account what is visible with CUDA only
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
        # HACK:
        equilibrationSteps = 100
        simulation.step(equilibrationSteps)

    # Add PLUMED bias
    with open(f'{output_dir}/{filename}_plumed.dat', 'r') as file:
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


