from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin
from sys import stdout
import fire

def main(
        filename=None
        ):
    if filename is None:
        raise ValueError('Filename is required')
    
    # 1: STRUCTURE
    # load the PDB system, can also load newer PDBx
    print('Loading...')
    pdb = PDBFile(f'{filename}_fixed.pdb')


    # 2: FORCE FIELD AND SYSTEM
    # load the forcefield that is specified in an xml file
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    # one can use more advanced AMBER forcefields with a more custom setup
    # or also CHARMM forcefields
    # e.g. we specifiy long-range interactions with PME (Particle Mesh Ewald)

    # invoke modeller to adjust the system
    modeller = Modeller(pdb.topology, pdb.positions)
    print('Adding hydrogens...')
    modeller.addHydrogens(forcefield, pH=7.0) # -> we fail here
    print('Adding solvent...')
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
    print('Platform used:', simulation.context.getPlatform().getName())
    # asign the positions
    simulation.context.setPositions(modeller.positions)
    # minimize the structure first (relaxation)
    print('Minimizing...')
    simulation.minimizeEnergy(maxIterations=1000)
    # # log the structure output as PDB
    # simulation.reporters.append(PDBReporter('output.pdb', 1000))
    # # log the energy and temperature every 1000 steps
    # simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
    #         # potentialEnergy=True, temperature=True))
    # # there's an automated way for checkpointing - CheckpointReporter

    # simulation.step(20000)

    print('Saving...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'{filename}_solvated.pdb', 'w'))

    print('Done')


if __name__ == '__main__':
    fire.Fire(main)