from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin
from sys import stdout
import fire


def main(
        filename=None,
        nsteps=5000
        ):
    
    if filename is None:
        raise ValueError('Filename is required')
    
    pdb = PDBFile(f'{filename}_solvated.pdb')

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    # log the structure output as PDB
    simulation.reporters.append(PDBReporter(f'{filename}_traj.pdb', 1000))
    # log the energy and temperature every 1000 steps
    simulation.reporters.append(
        StateDataReporter(
            # file=f'{filename}.out',
            file=stdout, 
            reportInterval=1000,
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
    # there's an automated way for checkpointing - CheckpointReporter


    import time

    t0 = time.time()
    simulation.step(nsteps)
    print('Total elapsed time:', time.time() - t0)
    print('Time per MD step', (time.time() - t0) / nsteps)



if __name__ == '__main__':
    fire.Fire(main)