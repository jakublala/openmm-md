from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin, nanosecond
from sys import stdout
import fire

# Define a function to check if a residue is a water molecule or an ion
def is_water_or_ion(residue):
    if residue.name == 'HOH' or residue.name in ['NA', 'CL']:
        return True
    return False



def main(
        filename=None,
        nsteps=None,
        mdtime=None, # give in ns
        ):
    
    time_step = 0.004*picoseconds

    if filename is None:
        raise ValueError('Filename is required')
    if nsteps is None and mdtime is None:
        raise ValueError('Number of steps or time is required')
    if nsteps is not None and mdtime is not None:
        if mdtime != time_step * nsteps:
            raise ValueError('Number of steps and time do not match')
    else:
        if mdtime is not None:
            nsteps = int(mdtime * nanosecond / time_step)
        elif nsteps is not None:
            mdtime = nsteps * time_step
        else:
            raise ValueError('Something went wrong...')

    

    
    pdb = PDBFile(f'tmp/{filename}_solvated.pdb')
    non_water_ion_atoms_indices = [atom.index for atom in pdb.topology.atoms() if not is_water_or_ion(atom.residue)]
    
    
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    system = forcefield.createSystem(
        topology=pdb.topology, 
        nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer,
        constraints=HBonds
        )
    print('Periodic boundary conditions:', system.usesPeriodicBoundaryConditions())

    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, time_step)

    
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': '0,1'}
    except:
        platform = Platform.getPlatformByName('OpenCL')
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    print('platform used:', simulation.context.getPlatform().getName())
    simulation.context.setPositions(pdb.positions)
    # log the structure output as PDB
    simulation.reporters.append(
        XTCReporter(
            f'tmp/{filename}.xtc', 
            # reportInterval=100,
            reportInterval=10000,
            enforcePeriodicBox=False,
            # atomSubset=non_water_ion_atoms_indices
            )
        )
    # log the energy and temperature every 1000 steps
    simulation.reporters.append(
        StateDataReporter(
            file=f'tmp/{filename}.out',
            # file=stdout,
            reportInterval=10000,
            # reportInterval=100,
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

    # move from tmp/ to output/
    import shutil
    shutil.move(f'tmp/{filename}.xtc', f'output/{filename}.xtc')
    shutil.move(f'tmp/{filename}.out', f'output/{filename}.out')
    shutil.move(f'tmp/{filename}_solvated.pdb', f'output/{filename}_solvated.pdb')



if __name__ == '__main__':
    fire.Fire(main)