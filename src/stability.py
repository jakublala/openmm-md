from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin, nanosecond
import fire
import time
from mdareporter import MDAReporter
from src.utils import get_checkpoint_interval, get_platform_and_properties

# Define a function to check if a residue is a water molecule or an ion
def is_water_or_ion(residue):
    if residue.name == 'HOH' or residue.name in ['NA', 'CL']:
        return True
    return False

def stability(
        filename=None,
        device_index='0',
        mdtime=None, # give in ns
        timestep=2, # in femtoseconds
        temperature=300,
        restart=False,
        logging_frequency=100,
        device="cuda",
        output_dir=None
        ):
    
    time_step = 0.001*timestep*picoseconds

    if output_dir is None:
        raise ValueError('Output directory is required')
    if filename is None:
        raise ValueError('Filename is required')
    
    nsteps = int(mdtime * nanosecond / time_step)
    
    if restart:
        pdb = PDBFile(f'output/{filename}/{filename}_solvated.pdb')
    else:
        pdb = PDBFile(f'{output_dir}/{filename}_solvated.pdb')
    
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    system = forcefield.createSystem(
        topology=pdb.topology, 
        nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer,
        constraints=HBonds
        )
    print('Periodic boundary conditions:', system.usesPeriodicBoundaryConditions())

    # TODO: run it at body temperature
    # have a bit larger friction coefficient initially in an equilibriation phase
    # then decrease it to a smaller value for production run
    integrator = LangevinMiddleIntegrator(temperature*kelvin, 1/picosecond, time_step)

    platform, properties = get_platform_and_properties(device, device_index)

    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    if restart:
        raise ValueError('Restart not implemented yet')
        simulation.loadCheckpoint(f'output/{filename}.chk')
    else:
        simulation.context.setPositions(pdb.positions)
    
    simulation.reporters.append(
        MDAReporter(
            f'{output_dir}/{filename}.dcd', 
            logging_frequency, 
            enforcePeriodicBox=False, 
            selection="protein"
            )
        )
    
    # log the energy and temperature every 1000 steps
    simulation.reporters.append(
        StateDataReporter(
            file=f'{output_dir}/{filename}.out',
            reportInterval=logging_frequency,
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

    simulation.reporters.append(
        CheckpointReporter(
            file=f'{output_dir}/{filename}.chk', 
            reportInterval=get_checkpoint_interval(timestep)
            )
        )

    t0 = time.time()
    simulation.step(nsteps)
    print('Total elapsed time:', time.time() - t0)
    print('Time per MD step', (time.time() - t0) / nsteps)

if __name__ == '__main__':
    fire.Fire(stability)