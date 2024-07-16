from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin
from sys import stdout
import fire

# Define a function to check if a residue is a water molecule or an ion
def is_water_or_ion(residue):
    if residue.name == 'HOH' or residue.name in ['NA', 'CL']:
        return True
    return False



def main(
        filename=None,
        nsteps=250000
        ):
    
    if filename is None:
        raise ValueError('Filename is required')
    
    pdb = PDBFile(f'{filename}_solvated.pdb')
    non_water_ion_atoms_indices = [atom.index for atom in pdb.topology.atoms() if not is_water_or_ion(atom.residue)]
    
    
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    system = forcefield.createSystem(
        topology=pdb.topology, 
        nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer,
        constraints=None
        )
    print('Periodic boundary conditions:', system.usesPeriodicBoundaryConditions())

    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    # log the structure output as PDB
    simulation.reporters.append(
        PDBReporter(
            f'{filename}_traj.pdb', 
            reportInterval=2,
            atomSubset=non_water_ion_atoms_indices
            )
        )
    # log the energy and temperature every 1000 steps
    simulation.reporters.append(
        StateDataReporter(
            file=f'{filename}.out',
            reportInterval=2,
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