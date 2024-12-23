from openmm import *
from openmm.app import *
from openmm.unit import *
from openmmplumed import PlumedForce
from openmm.unit import nanometers, picoseconds, kelvin, nanoseconds
from mdareporter import MDAReporter
from openmm.app import PDBxFile
import logging
import os
logger = logging.getLogger(__name__)

from src.utils import get_checkpoint_interval, get_platform_and_properties
from src.plumed.io import create_plumed_input
from src.utils import save_equilibrated_state

def run_plumed(
        filename, 
        mdtime, 
        timestep, 
        temperature=300, 
        restart_checkpoint=None, 
        device='cuda', 
        device_precision='mixed',
        device_index='0', 
        output_dir=None, 
        logging_frequency=None,
        plumed_config=None,
        plumed_mode=None,
        equilibrate_only=False,
        generate_plumed_input=True,
        ):
    """Run an OpenMM molecular dynamics simulation with Metadynamics or OPES (On-the-fly Probability Enhanced Sampling).
    
    Parameters
    ----------
    filename : str
        Base name for input/output files. The function expects '{filename}_solvated.cif' as input
        and will generate various output files using this base name.
    mdtime : float
        Total simulation time in nanoseconds.
    timestep : float
        Integration timestep multiplier. The actual timestep will be 0.001 * timestep picoseconds.
    device_index : str, default='0'
        GPU device index to use for CUDA calculations.
    temperature : float, default=300
        Simulation temperature in Kelvin.
    restart_checkpoint : str, optional
        Path to a checkpoint file to restart simulation from.
    device : str, default='cuda'
        Computation device to use. Options are 'cuda', 'cpu', or 'opencl'.
    output_dir : str, required
        Directory path for all input and output files.
    logging_frequency : float, required
        Frequency (in picoseconds) for trajectory and state data reporting.
    plumed_config : dict, required
        PLUMED config.
    plumed_mode : str, required
        Chain mode

    Returns
    -------
    None
        The function saves various output files to the specified output directory:
        - {filename}.dcd: Trajectory file
        - {filename}.out: Simulation statistics and progress
        - {filename}.chk: Checkpoint files
        
    Notes
    -----
    The simulation uses the AMBER14 force field with TIP3P water model and includes:
    - PME (Particle Mesh Ewald) for long-range interactions
    - 1.0 nm nonbonded cutoff
    - Hydrogen bond constraints
    - Rigid water molecules
    - Langevin middle integrator for temperature control
    
    The simulation includes an equilibration phase (1 ns) if not starting from a checkpoint,
    and applies PLUMED bias according to a script specified in '{filename}_plumed.dat'.

    Raises
    ------
    ValueError
        If output_dir is not specified or if an invalid device is selected.
    """

    if output_dir is None:
        raise ValueError('Output directory is required')
    
    if os.path.exists(f'{output_dir}/{filename}_equilibrated.cif'):
        logger.info(f'Equilibrated state found at {output_dir}/{filename}_equilibrated.cif')
        logger.info('Skipping equilibration...')
        pdb = PDBxFile(f'{output_dir}/{filename}_equilibrated.cif')
        equilibrated = True
    else:
        pdb = PDBxFile(f'{output_dir}/{filename}_solvated.cif')
        equilibrated = False

    
    if not plumed_config['trajectory_logging']:
        plumed_config['trajectory_logging'] = True


    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # System Configuration
    nonbondedMethod = PME
    nonbondedCutoff = 1.0*nanometers  # This cutoff and the following for tolerance are standard values
    ewaldErrorTolerance = 10**(-4)
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 10**(-4) # quite a default value for accuracy

    # Integration Options
    dt = 0.001 * picoseconds * timestep
    temperature = temperature * kelvin
    friction = 1.0/picoseconds

    # Simulation Options
    steps = int(mdtime * nanoseconds / dt)
    equilibrationSteps = int(1 * nanoseconds / dt)

    traj_interval = int(logging_frequency * picoseconds / dt)

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
    checkpointReporter = CheckpointReporter(
        f'{output_dir}/{filename}.chk', 
        get_checkpoint_interval(timestep)
        )


    # Prepare the Simulation
    logger.info('Building system...')
    topology = pdb.topology
    positions = pdb.positions
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
    
    platform, properties = get_platform_and_properties(device, device_index, device_precision)

    simulation = Simulation(topology, system, integrator, platform, properties)

    # get the periodic box vectors and log them
    box_vectors = simulation.topology.getPeriodicBoxVectors()

    # Print the box vectors
    logger.info("Box Vectors (in nanometers):")
    logger.info(f"x: {box_vectors[0] / nanometers}")
    logger.info(f"y: {box_vectors[1] / nanometers}")
    logger.info(f"z: {box_vectors[2] / nanometers}")
    
    if restart_checkpoint: 
        logger.info(f'Restarting the OpenMM simulation from a checkpoint at {restart_checkpoint}')
        simulation.loadCheckpoint(restart_checkpoint)
        # no equilibration for system from checkpoint
    elif equilibrated:
        if equilibrate_only:
            raise ValueError("Equilibrate_only is True, but equilibrated state found. You don't need to equilibrate again.")
        simulation.context.setPositions(positions)
        simulation.context.setVelocitiesToTemperature(temperature)
    else:
        simulation.context.setPositions(positions)

        # Equilibrate
        logger.info('Equilibrating...')
        simulation.context.setVelocitiesToTemperature(temperature)
        # TODO: make it random, by having the number of steps be also somewhat number
        simulation.step(equilibrationSteps)

        save_equilibrated_state(
            simulation=simulation,
            output_dir=output_dir,
            filename=filename
            )
        
        if equilibrate_only:
            return
        
    if plumed_config is None:
        raise ValueError('PLUMED config is required')
    
    if generate_plumed_input:
        create_plumed_input(
            filename=filename, 
            output_dir=output_dir,
            config=plumed_config,
            mode=plumed_mode
            )
    else:
        if not os.path.exists(f'{output_dir}/{filename}_plumed.dat'):
            raise ValueError(f'PLUMED input file {output_dir}/{filename}_plumed.dat not found')
        else:
            logger.info(f'PLUMED input file {output_dir}/{filename}_plumed.dat found, skipping generation. Possibly a manual overwrite.')

    # Add PLUMED bias
    with open(f'{output_dir}/{filename}_plumed.dat', 'r') as file:
        script = file.read()

    plumed_force = PlumedForce(script)
    system.addForce(plumed_force)

    # Before creating the simulation, assign unique force groups
    # This might not be necessary, but it's a good idea to ensure visibility of force groups
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)
    
    # Reinitialize the simulation with the updated system
    simulation.system = system
    simulation.context.reinitialize(preserveState=True)

    # Simulate
    logger.info('Simulating...')

    if plumed_config['trajectory_logging']:
        simulation.reporters.append(trajReporter)
    simulation.reporters.append(dataReporter)
    simulation.reporters.append(checkpointReporter)
    simulation.currentStep = 0
    simulation.step(steps)    

    # TODO: DO ADAPTIVE CONVERGENCE HERE!!!!!
    # i.e. extend the simulation if convergence is not reached
