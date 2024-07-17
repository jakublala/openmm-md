from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin, kilojoules_per_mole
from sys import stdout
import fire
import logging


logger = logging.getLogger(__name__)

def main(
        filename=None
        ):
    if filename is None:
        raise ValueError('Filename is required')
    
    # 1: STRUCTURE
    # load the PDB system, can also load newer PDBx
    print('Loading...', flush=True)
    pdb = PDBFile(f'{filename}_fixed.pdb')


    # 2: FORCE FIELD AND SYSTEM
    # load the forcefield that is specified in an xml file
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    # one can use more advanced AMBER forcefields with a more custom setup
    # or also CHARMM forcefields
    # e.g. we specifiy long-range interactions with PME (Particle Mesh Ewald)

    # invoke modeller to adjust the system
    modeller = Modeller(pdb.topology, pdb.positions)
    print('Adding hydrogens...', flush=True)
    modeller.addHydrogens(forcefield, pH=7.0) # -> we fail here
    print('Adding solvent...', flush=True)
    modeller.addSolvent(forcefield, padding=1*nanometer)

    # modeller could also add a membrane or ions (ionic strength of the solvent)
    # print('Adding extra (virtual) particles???')
    # modeller.addExtraParticles()




    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)

    # 3: INTEGRATOR (how do we run the dynamics, what ensemble do we run in)
    # integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    integrator = VerletIntegrator(0.001*picoseconds)

    # 4: SIMULATION
    # there's 4 different platforms tu use - Reference, CPU, CUDA, and OpenCL
    # platform = Platform.getPlatformByName('CUDA')
    # then assign it to Simulation object

    simulation = Simulation(modeller.topology, system, integrator)
    print('Platform used:', simulation.context.getPlatform().getName(), flush=True)
    # asign the positions
    simulation.context.setPositions(modeller.positions)
    # minimize the structure first (relaxation)
    import time
    class Reporter(MinimizationReporter):
        interval = 10 # report interval
        iterations = []
        energies = [] # array to record progress
        times = [time.time()] # array to record progress

        def report(self, iteration, x, grad, args):
            # print current system energy to screen 
            if iteration % self.interval == 0:
                print(iteration, args['system energy'], time.time()-self.times[-1], flush=True)

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

    print("Minimizing energy...", flush=True)
    simulation.minimizeEnergy(
        tolerance=10*kilojoules_per_mole/nanometer,
        maxIterations=100,
        reporter=Reporter(),
    )

    print('Saving...', flush=True)
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'{filename}_solvated.pdb', 'w'))

    # import matplotlib.pyplot as plt
    # # plot energies
    # for i, (iteration, energy) in enumerate(zip(Reporter.iterations, Reporter.energies)):
    #     plt.plot(iteration, energy, 'x', alpha=(i+1)/len(Reporter.energies))
    # plt.xlabel('Iteration')
    # plt.ylabel('Energy')
    # plt.savefig(f'{filename}_minimization.png')

    print('Done', flush=True)


if __name__ == '__main__':
    fire.Fire(main)