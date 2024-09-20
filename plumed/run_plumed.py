from openmm import *
from openmm.app import *
from openmm.unit import *
from openmm import unit

import random
import numpy as np
random.seed(42)
np.random.seed(42)

pdf = PDBFile('equilibrated.pdb')

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# System Configuration
nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers  #This cutoff and the following for tolerance are standard values
ewaldErrorTolerance = 10**(-4)
constraints = HBonds
rigidWater = True
constraintTolerance = 10**(-4) #Quite a default value for accuracy

# Integration Options
fs = 0.001 * picoseconds
dt = 4 * fs
temperature = 300*kelvin
friction = 1.0/picosecond

# Simulation Options
steps = 10**7
equilibrationSteps = 10**5
# platform = Platform.getPlatformByName('OpenCL')
platform = Platform.getPlatformByName('CUDA')


# TODO: employ this later on with 100 ps trajectory
traj_interval = int(100 * picoseconds / dt)

dcdReporter = DCDReporter(
	'traj.dcd', 
	10000,
	enforcePeriodicBox=False
	)

dataReporter = StateDataReporter(
            file='logger',
            # reportInterval=10000,
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
            totalSteps=steps,
			separator='\t'
        )
checkpointReporter = CheckpointReporter('checkpoint', 10000)


# Prepare the Simulation
print('Building system...')
topology = pdf.topology
positions = pdf.positions
system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)


# Add PLUMED bias
with open('plumed.dat', 'r') as file:
    script = file.read()

from openmmplumed import PlumedForce
plumed_force = PlumedForce(script)
system.addForce(plumed_force)

integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)

# Simulate
print('Simulating...')

simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
simulation.step(steps)
