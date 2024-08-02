from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin, kilojoules_per_mole, nanosecond, bar, amu
from sys import stdout
import fire
import logging
import numpy as np


from utils import get_full_reporter

logger = logging.getLogger(__name__)

def minimize(
        filename=None,
        max_iterations=None,
        device_index='0',
        constraints=None,
        ):
    """
    Uses the OpenMM library (i.e. LBFGS) to minimize the energy of a given PDB file.
    
    """

    logging.info('Minimizing energy...')

    if filename is None:
        raise ValueError('Filename is required')
    if max_iterations is None:
        raise ValueError('Max iterations is required')


    # 1: STRUCTURE
    # load the PDB system, can also load newer PDBx
    logger.info('Loading PDB file...')
    pdb = PDBFile(f'tmp/{filename}_fixed.pdb')


    # 2: FORCE FIELD AND SYSTEM
    # load the forcefield that is specified in an xml file
    forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml')
    # one can use more advanced AMBER forcefields with a more custom setup
    # or also CHARMM forcefields
    # e.g. we specifiy long-range interactions with PME (Particle Mesh Ewald)

    # invoke modeller to adjust the system
    modeller = Modeller(pdb.topology, pdb.positions)
    logger.info('Adding hydrogens...')
    modeller.addHydrogens(forcefield, pH=7.0) # -> we fail here
    logger.info('Adding solvent...')
    modeller.addSolvent(forcefield, padding=1*nanometer)

    # modeller could also add a membrane or ions (ionic strength of the solvent)

    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=constraints)

    # 3: INTEGRATOR (how do we run the dynamics, what ensemble do we run in)
    # integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    integrator = VerletIntegrator(0.001*picoseconds)

    # 4: SIMULATION
    # there's 4 different platforms tu use - Reference, CPU, CUDA, and OpenCL
    # platform = Platform.getPlatformByName('CUDA')
    # then assign it to Simulation object
    properties = None
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': device_index}
    except:
        platform = Platform.getPlatformByName('OpenCL')


    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    logger.info(f'Platform used: {simulation.context.getPlatform().getName()}')
    # asign the positions
    simulation.context.setPositions(modeller.positions)
    # minimize the structure first (relaxation)
    import time
    class Reporter(MinimizationReporter):
        interval = 100 # report interval
        iterations = []
        energies = [] # array to record progress
        times = [time.time()] # array to record progress

        def report(self, iteration, x, grad, args):
            # save energy at each iteration to an array we can use later
            if iteration == 0: # create new list
                self.energies.append([args['system energy']])
                self.iterations.append([iteration])
            else:
                self.energies[-1].append(args['system energy'])
                self.iterations[-1].append(iteration)

            self.times.append(time.time())

            # The report method must return a bool specifying if minimization should be stopped. 
            # You can use this functionality for early termination.
            return False

    logger.info("Minimizing energy...")
    simulation.minimizeEnergy(
        tolerance=5*kilojoules_per_mole/nanometer,
        maxIterations=max_iterations,
        reporter=Reporter(),
    )

    logger.info('Saving...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'tmp/{filename}_solvated.pdb', 'w'))

    import matplotlib.pyplot as plt
    # plot energies
    for i, (iteration, energy) in enumerate(zip(Reporter.iterations, Reporter.energies)):
        plt.plot(iteration, energy, 'b-', alpha=(i+1)/len(Reporter.energies), label=f'Iteration {i+1}')
    plt.xlabel('LBFGS Step')
    plt.ylabel('Energy / kJ/mole')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'tmp/{filename}_minimization.png', dpi=300)
    plt.clf()

    logger.info('Done')



def relax_md_npt(
        filename=None,
        mdtime = None,
        device_index='0',
        constraints=None,
        fix=None,
        annealing=True,
):
    """
    Runs MD in NPT to relax the system.
    """
    logger.info('Equilibrating with NPT (relax_md_npt)...')

    if annealing:
        current_T = 100
    else:
        current_T = 300 # TODO: check that this system works at this temp (or whether I could be at body temperature)
    final_T = 300
    timestep = 0.004*picoseconds / np.sqrt(12)
    if mdtime is None:
        mdtime = 0.01*nanosecond # ns
        logger.warning(f'mdtime is not set, using default value of {mdtime}')
    nsteps = int(mdtime / timestep)
    

    if filename is None:
        raise ValueError('Filename is required')

    # 1: STRUCTURE
    pdb = PDBFile(f'tmp/{filename}_solvated.pdb')

    # 2: FORCE FIELD AND SYSTEM
    forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=constraints)
    
    if fix is not None:
        if fix == 'protein':
            logger.info('Fixing protein atoms...')
            for atom in pdb.topology.atoms():
                if atom.residue.name != 'HOH':
                    system.setParticleMass(atom.index, 0*amu)
        elif fix == 'water':
            logger.info('Fixing water atoms...')
            for atom in pdb.topology.atoms():
                if atom.residue.name == 'HOH':
                    system.setParticleMass(atom.index, 0*amu)
        else:
            raise ValueError('Invalid value for fix')
    else:
        logger.info('Not fixing any atoms...')
        fix = 'all'


    logger.info(f'Periodic boundary conditions: {system.usesPeriodicBoundaryConditions()}')

    # 3: INTEGRATOR (NPT ensemble)
    integrator = LangevinIntegrator(current_T*kelvin, 1/picosecond, timestep)
    system.addForce(MonteCarloBarostat(1*bar, current_T*kelvin))

    # 4: SIMULATION
    properties = None
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': device_index}
    except:
        platform = Platform.getPlatformByName('OpenCL')

    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(current_T*kelvin)
    simulation.reporters.append(
        get_full_reporter(
            f"{filename}_npt_{fix}", 
            log_freq=1000, 
            nsteps=nsteps
            )
        )

    # 5: EQUILIBRATE
    # annealing from 100 K to 300 K
    if annealing:
        adjust_freq = 1000
        logger.info('Annealing...')
        logger.info(f'{nsteps=}, {adjust_freq=}, {timestep=}, {current_T=}, {final_T=}')
        for i in range(adjust_freq):
            current_T = (current_T + i*(final_T - current_T)/adjust_freq) * kelvin
            integrator.setTemperature(current_T)
            simulation.context.setParameter(MonteCarloBarostat.Temperature(), current_T)
            simulation.step(nsteps//adjust_freq)
    else:
        logger.info('Equilibriating...')
        simulation.step(nsteps)


    # save the final structure
    logger.info('Saving...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'tmp/{filename}_solvated.pdb', 'w'))




if __name__ == '__main__':
    fire.Fire(minimize)