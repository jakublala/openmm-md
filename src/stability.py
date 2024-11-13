from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, picoseconds, kelvin, nanosecond
from sys import stdout
import fire
import time
from mdareporter import MDAReporter

# Define a function to check if a residue is a water molecule or an ion
def is_water_or_ion(residue):
    if residue.name == 'HOH' or residue.name in ['NA', 'CL']:
        return True
    return False

def stability(
        filename=None,
        nsteps=None,
        mdtime=None, # give in ns
        restart=False,
        device_index='0',
        log_freq=10000,
        timestep=4, # in femtoseconds
        device="cuda",
        output_dir=None
        ):
    
    time_step = 0.001*timestep*picoseconds

    if output_dir is None:
        raise ValueError('Output directory is required')
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

    log_freq = int(100 * picoseconds / time_step)

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
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, time_step)

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
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    print('platform used:', simulation.context.getPlatform().getName())
    if restart:
        raise ValueError('Restart not implemented yet')
        simulation.loadCheckpoint(f'output/{filename}.chk')
    else:
        simulation.context.setPositions(pdb.positions)
    
    simulation.reporters.append(
        MDAReporter(
            f'{output_dir}/{filename}.dcd', 
            log_freq, 
            enforcePeriodicBox=False, 
            selection="protein"
            )
        )
    
    # log the energy and temperature every 1000 steps
    simulation.reporters.append(
        StateDataReporter(
            file=f'{output_dir}/{filename}.out',
            reportInterval=log_freq,
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
            reportInterval=log_freq*1000000000
            )
        )

    t0 = time.time()
    simulation.step(nsteps)
    print('Total elapsed time:', time.time() - t0)
    print('Time per MD step', (time.time() - t0) / nsteps)

if __name__ == '__main__':
    fire.Fire(stability)