from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin, kilojoules_per_mole, nanosecond, bar, amu
from sys import stdout
import fire
import logging
import numpy as np


from src.utils import get_full_reporter

logger = logging.getLogger(__name__)

def minimize(
        filename=None,
        max_iterations=None,
        device_index='0',
        constraints=None,
        device="cuda",
        ):
    """
    Uses the OpenMM library (i.e. LBFGS) to minimize the energy of a given PDB file.
    
    """

    if filename is None:
        raise ValueError('Filename is required')
    if max_iterations is None:
        raise ValueError('Max iterations is required')


    # 1: STRUCTURE
    # load the PDB system, can also load newer PDBx
    logger.info('Loading PDB file...')
    pdb_file = PDBFile(f'tmp/{filename}/{filename}_fixed.pdb')


    # 2: FORCE FIELD AND SYSTEM
    # load the forcefield that is specified in an xml file
    forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml')
    # one can use more advanced AMBER forcefields with a more custom setup
    # or also CHARMM forcefields
    # e.g. we specifiy long-range interactions with PME (Particle Mesh Ewald)

    # invoke modeller to adjust the system
    modeller = Modeller(pdb_file.topology, pdb_file.positions)
    logger.info('Adding hydrogens...')
    modeller.addHydrogens(forcefield, pH=7.0)
    logger.info('Adding solvent...')
    modeller.addSolvent(forcefield, padding=1*nanometer)

    # modeller could also add a membrane or ions (ionic strength of the solvent)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=constraints)

    # 3: INTEGRATOR (how do we run the dynamics, what ensemble do we run in)
    # integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

    # 4: SIMULATION
    # there's 4 different platforms tu use - Reference, CPU, CUDA, and OpenCL
    # platform = Platform.getPlatformByName('CUDA')
    # then assign it to Simulation object
    properties = None
    if device == "cuda":
        try:
            platform = Platform.getPlatformByName('CUDA')
            properties = {'DeviceIndex': device_index}
        except:
            print('CUDA not available, using OpenCL')
            platform = Platform.getPlatformByName('OpenCL')
    elif device == "cpu":
        platform = Platform.getPlatformByName('OpenCL')
    else:
        raise ValueError('Invalid device')

    logger.info(f'Platform used: {platform.getName()}')
    
    import time
    class Reporter(MinimizationReporter):
        def __init__(self):
            super(Reporter, self).__init__()
            self.interval = 100 # report interval
            self.iterations = []
            self.energies = [] # array to record progress
            self.times = [time.time()] # array to record progress

        def report(self, iteration, x, grad, args):
            # if iteration % self.interval == 0:
            #     logger.info(f'{self.__class__.__name__}: Iteration {iteration}, Energy {args["system energy"]}, Grad {np.linalg.norm(grad)}')

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
        
    import matplotlib.pyplot as plt

    # Create separate Reporter instances for each minimization
    reporter_hyrogens = Reporter()
    reporter_water = Reporter()
    reporter_protein = Reporter()
    reporter_whole = Reporter()

    protein_ids = [atom.index for atom in modeller.topology.atoms() if atom.residue.name != 'HOH']
    water_ids = [atom.index for atom in modeller.topology.atoms() if atom.residue.name == 'HOH']
    non_hydrogen_ids = [atom.index for atom in modeller.topology.atoms() if atom.element.symbol != 'H']
    
    # save the topology now into a PDB
    PDBFile.writeFile(modeller.topology, modeller.positions, open(f'tmp/{filename}/{filename}_solvated.pdb', 'w'))

    # 0: first first minimize just the hydrogens
    non_hydrogen_restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    system.addForce(non_hydrogen_restraint)
    non_hydrogen_restraint.addGlobalParameter('k', 1e10*kilojoules_per_mole/nanometer)
    non_hydrogen_restraint.addPerParticleParameter('x0')
    non_hydrogen_restraint.addPerParticleParameter('y0')
    non_hydrogen_restraint.addPerParticleParameter('z0')
    for id in non_hydrogen_ids:
        non_hydrogen_restraint.addParticle(id, modeller.positions[id])
    logger.info("Minimizing energy (LBFGS) with fixed non-hydrogen atoms...")

    integrator = VerletIntegrator(0.001*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    original_pos = simulation.context.getState(getPositions=True).getPositions()
    original_pos = np.array(original_pos.value_in_unit(nanometer))
    simulation.minimizeEnergy(
        tolerance=10*kilojoules_per_mole/nanometer,
        maxIterations=max_iterations,
        reporter=reporter_hyrogens,
    )
    system.removeForce(system.getNumForces()-1)
    
    new_pos = simulation.context.getState(getPositions=True).getPositions()
    new_pos = np.array(new_pos.value_in_unit(nanometer))
    
    # assert that protein and water positions are roughly the same within some tolerance (non-hydrogens)
    hydrogen_ids = [atom.index for atom in modeller.topology.atoms() if atom.element.symbol == 'H']

    hydrogen_dists = np.linalg.norm(new_pos[hydrogen_ids] - original_pos[hydrogen_ids], axis=1)
    non_hydrogen_dists = np.linalg.norm(new_pos[non_hydrogen_ids] - original_pos[non_hydrogen_ids], axis=1)

    # Determine the overall range
    min_val = min(np.min(hydrogen_dists), np.min(non_hydrogen_dists))
    max_val = max(np.max(hydrogen_dists), np.max(non_hydrogen_dists))

    # Create bins based on this range
    bins = np.linspace(min_val, max_val, 101)  # 101 edges to create 100 bins

    # Plot the histograms
    plt.hist(hydrogen_dists, bins=bins, alpha=0.5, label='hydrogens', density=True)
    plt.hist(non_hydrogen_dists, bins=bins, alpha=0.5, label='non-hydrogens', density=True)
    plt.xlabel('Distance (nm)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title('Distance moved by atoms during water relaxation')
    plt.savefig(f'tmp/{filename}/{filename}_nonhydrogens_fixed.png', dpi=300)
    plt.clf()

    assert np.all(non_hydrogen_dists < 0.1), 'Protein atoms moved during water relaxation (moved more than 0.1 nm)'
    del simulation, integrator, new_pos, original_pos

    # A: first minimize just the water
    protein_restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    system.addForce(protein_restraint)
    protein_restraint.addGlobalParameter('k', 1e10*kilojoules_per_mole/nanometer)
    protein_restraint.addPerParticleParameter('x0')
    protein_restraint.addPerParticleParameter('y0')
    protein_restraint.addPerParticleParameter('z0')
    for id in protein_ids:
        protein_restraint.addParticle(id, modeller.positions[id])
    logger.info("Minimizing energy (LBFGS) with fixed protein atoms...")

    integrator = VerletIntegrator(0.001*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    original_pos = simulation.context.getState(getPositions=True).getPositions()
    original_pos = np.array(original_pos.value_in_unit(nanometer))
    simulation.minimizeEnergy(
        tolerance=10*kilojoules_per_mole/nanometer,
        maxIterations=max_iterations,
        reporter=reporter_water,
    )
    system.removeForce(system.getNumForces()-1)
    
    new_pos = simulation.context.getState(getPositions=True).getPositions()
    new_pos = np.array(new_pos.value_in_unit(nanometer))
    
    # assert that protein positions are roughly the same within some tolerance
    
    water_dists = np.linalg.norm(new_pos[water_ids] - original_pos[water_ids], axis=1)
    protein_dists = np.linalg.norm(new_pos[protein_ids] - original_pos[protein_ids], axis=1)

    # Determine the overall range
    min_val = min(np.min(water_dists), np.min(protein_dists))
    max_val = max(np.max(water_dists), np.max(protein_dists))

    # Create bins based on this range
    bins = np.linspace(min_val, max_val, 101)  # 101 edges to create 100 bins

    # Plot the histograms
    plt.hist(water_dists, bins=bins, alpha=0.5, label='water', density=True)
    plt.hist(protein_dists, bins=bins, alpha=0.5, label='protein', density=True)
    plt.xlabel('Distance (nm)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title('Distance moved by atoms during water relaxation')
    plt.savefig(f'tmp/{filename}/{filename}_protein_fixed.png', dpi=300)
    plt.clf()

    assert np.all(protein_dists < 0.1), 'Protein atoms moved during water relaxation (moved more than 0.1 nm)'
    del simulation, integrator, new_pos, original_pos


    # # B: then minimize just the protein
    # water_restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    # system.addForce(water_restraint)
    # water_restraint.addGlobalParameter('k', 1e10*kilojoules_per_mole/nanometer)
    # water_restraint.addPerParticleParameter('x0')
    # water_restraint.addPerParticleParameter('y0')
    # water_restraint.addPerParticleParameter('z0')
    # for id in water_ids:
    #     water_restraint.addParticle(id, modeller.positions[id])
    # logger.info("Minimizing energy (LBFGS) with fixed water atoms...")
    # integrator = VerletIntegrator(0.001*picoseconds)
    # simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    # simulation.context.setPositions(modeller.positions)
    # original_pos = simulation.context.getState(getPositions=True).getPositions()
    # original_pos = np.array(original_pos.value_in_unit(nanometer))
    # simulation.minimizeEnergy(
    #     tolerance=5*kilojoules_per_mole/nanometer,
    #     maxIterations=max_iterations,
    #     reporter=reporter_protein,
    # )
    # simulation.system.removeForce(system.getNumForces()-1)

    # new_pos = simulation.context.getState(getPositions=True).getPositions()
    # new_pos = np.array(new_pos.value_in_unit(nanometer))
    
    # # assert that protein positions are roughly the same within some tolerance
    # water_dists = np.linalg.norm(new_pos[water_ids] - original_pos[water_ids], axis=1)
    # protein_dists = np.linalg.norm(new_pos[protein_ids] - original_pos[protein_ids], axis=1)

    # # Determine the overall range
    # min_val = min(np.min(water_dists), np.min(protein_dists))
    # max_val = max(np.max(water_dists), np.max(protein_dists))

    # # Create bins based on this range
    # bins = np.linspace(min_val, max_val, 101)  # 101 edges to create 100 bins

    # # Plot the histograms
    # plt.hist(water_dists, bins=bins, alpha=0.5, label='water', density=True)
    # plt.hist(protein_dists, bins=bins, alpha=0.5, label='protein', density=True)
    # plt.xlabel('Distance (nm)')
    # plt.ylabel('Frequency')
    # plt.legend()
    # plt.title('Distance moved by atoms during protein relaxation')
    # plt.savefig(f'tmp/{filename}_water_fixed.png', dpi=300)
    # plt.clf()

    # assert np.all(water_dists < 0.001), 'Water atoms moved during protein relaxation (moved more than 0.001 nm)'
    # del simulation, integrator, new_pos, original_pos

    # C: then minimize the whole system
    logger.info("Minimizing energy (LBFGS) on the whole system (protein + water)...")
    integrator = VerletIntegrator(0.001*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    original_pos = simulation.context.getState(getPositions=True).getPositions()
    original_pos = np.array(original_pos.value_in_unit(nanometer))
    simulation.minimizeEnergy(
        tolerance=10*kilojoules_per_mole/nanometer,
        maxIterations=max_iterations, # TODO: remove this and see what happens
        reporter=reporter_whole,
    )

    new_pos = simulation.context.getState(getPositions=True).getPositions()
    new_pos = np.array(new_pos.value_in_unit(nanometer))
    
    water_dists = np.linalg.norm(new_pos[water_ids] - original_pos[water_ids], axis=1)
    protein_dists = np.linalg.norm(new_pos[protein_ids] - original_pos[protein_ids], axis=1)

    plt.hist(water_dists, bins=100, alpha=0.5, label='water', density=True)
    plt.hist(protein_dists, bins=100, alpha=0.5, label='protein', density=True)
    plt.xlabel('Distance (nm)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title('Distance moved by protein atoms during both (protein + water) relaxation')
    plt.savefig(f'tmp/{filename}/{filename}_both.png', dpi=300)
    plt.clf()

    logger.info('Saving...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'tmp/{filename}/{filename}_solvated.pdb', 'w'))

    # Plot energies
    plt.figure(figsize=(10, 6))

    colors = ['yellow', 'blue', 'green', 'red']
    labels = ['Hydrogens Relax', 'Water Relax', 'Protein Relax', 'Both Relax']

    for i, reporter in enumerate([reporter_hyrogens, reporter_water, reporter_protein, reporter_whole]):
        if len(reporter.iterations) == 0:
            continue
        for j, (iteration, energy) in enumerate(zip(reporter.iterations, reporter.energies)):
            alpha = (j+1) / len(reporter.energies)
            plt.plot(iteration, energy, color=colors[i], alpha=alpha)
        
        # Add solid line and label for the final iteration
        plt.plot(reporter.iterations[-1], reporter.energies[-1], color=colors[i], label=labels[i])

    plt.xlabel('LBFGS Step')
    plt.ylabel('Energy (kJ/mol)')
    plt.legend()
    plt.title('Energy Minimization Progress')
    plt.tight_layout()
    plt.savefig(f'tmp/{filename}/{filename}_minimization.png', dpi=300)
    plt.close()

    logger.info('Done')


def relax_md_npt(
        filename=None,
        mdtime = None,
        device_index='0',
        constraints=None,
        fix=None,
        annealing=True,
        timestep=4 # in femtoseconds
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
    if constraints is None:
        timestep = 0.001*timestep*picoseconds / np.sqrt(12)
    else: # assume HBonds constraint included, but assert
        assert constraints == HBonds, 'Only HBonds constraint is supported'
        timestep = 0.001*timestep*picoseconds
    if mdtime is None:
        mdtime = 0.01*nanosecond # ns
        logger.warning(f'mdtime is not set, using default value of {mdtime}')
    else:
        mdtime = mdtime * nanosecond
    nsteps = int(mdtime / timestep)
    

    if filename is None:
        raise ValueError('Filename is required')

    # 1: STRUCTURE
    pdb_file = PDBFile(f'tmp/{filename}_solvated.pdb')

    # 2: FORCE FIELD AND SYSTEM
    forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml')

    system = forcefield.createSystem(pdb_file.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=constraints)
    
    if fix is not None:
        if fix == 'protein':
            logger.info('Fixing protein atoms...')
            for atom in pdb_file.topology.atoms():
                if atom.residue.name != 'HOH':
                    system.setParticleMass(atom.index, 0*amu)
        elif fix == 'water':
            logger.info('Fixing water atoms...')
            for atom in pdb_file.topology.atoms():
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

    simulation = Simulation(pdb_file.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb_file.positions)
    simulation.context.setVelocitiesToTemperature(current_T*kelvin)

    if annealing:
        output_file = f'tmp/{filename}_npt_annealing_{fix}.out'
    else:
        output_file = f'tmp/{filename}_npt_equilibriating_{fix}.out'

    simulation.reporters.append(
        get_full_reporter(
            f"{filename}_npt_{fix}_{str(constraints)}", 
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
            # current_T here is actually the initial T
            T = (current_T + i*(final_T - current_T)/adjust_freq) * kelvin
            integrator.setTemperature(T)
            simulation.context.setParameter(MonteCarloBarostat.Temperature(), T)
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