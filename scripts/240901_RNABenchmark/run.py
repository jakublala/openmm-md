# This script was generated by OpenMM-Setup on 2024-08-18.

from openmm import *
from openmm.app import *
from openmm.unit import *
from openmm import unit
from pdbfixer import PDBFixer

from openmm import MonteCarloBarostat
import openmm as mm

from openmm.app import ForceField, Modeller, PDBFile
from openmm import unit

# Load the PDB file
fixer = PDBFixer(filename='6ol3.pdb')

# Load the force field
forcefield = ForceField('amber14/RNA.OL3.xml', 'amber14/tip3pfb.xml')

#STEP 1: CLEAN UP THE PDB
# Add missing residues (completes existing residues with missing atoms, especially terminal residues...)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
## Add missing hydrogens
fixer.addMissingHydrogens()


#STEP 2: SOLVATE THE STRUCTURE

# Create a Modeller object

# Jakub: open fixed pdb file

# modeller = Modeller( fixer.topology, fixer.positions )
# pdb_file = PDBFile('tmp/S2_Best_AB_fixed.pdb')
# modeller = Modeller(pdb_file.topology, pdb_file.positions)

Save the solvated structure
with open('tmp/fixed-unsolvated.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
# Create a Modeller object
modeller = Modeller(pdb.topology, pdb.positions)

# Define the box size and add water
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometers)

# Get the box vectors (i.e., box size)
box_vectors = modeller.topology.getPeriodicBoxVectors()

# Print the box vectors
print("Box Vectors (in nanometers):")
print("x:", box_vectors[0] / unit.nanometers)
print("y:", box_vectors[1] / unit.nanometers)
print("z:", box_vectors[2] / unit.nanometers)

# Save the solvated structure
with open('tmp/solvated.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)


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
steps = 10**6
equilibrationSteps = 10**5
platform = Platform.getPlatformByName('CUDA')
dcdReporter = DCDReporter('traj.out', 10)
dataReporter = StateDataReporter('logger',10, totalSteps=steps,
    step=True, time=True, speed=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter('checkpoint', 10000)

# Prepare the Simulation
print('Building system...')
topology = modeller.topology
positions = modeller.positions
system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)

# NOT USED ANYMORE SINCE THIS IS SET BY THE ADDSOLVENT PART
## Define the box dimensions
#box_length = 6. * nanometers
## Define box vectors
#box_vectors = [
#    mm.Vec3(box_length, 0, 0),
#    mm.Vec3(0, box_length, 0),
#    mm.Vec3(0, 0, box_length)
#]
## Set periodic box vectors
#system.setDefaultPeriodicBoxVectors(*box_vectors)


integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platform, {'DeviceIndex': '1'})
simulation.context.setPositions(positions)
# Optionally, add a barostat to the system
barostat = MonteCarloBarostat(1 * unit.atmosphere, 300 * unit.kelvin)
system.addForce(barostat)

## Minimize and Equilibrate. First, minimize ONLY positions of H atoms
## Create a dictionary to map atom indices to their types
atom_types = {}
for atom in topology.atoms():
    atom_types[atom.index] = atom.element.symbol  # Use element symbol or type
#
## Identify hydrogen atoms
hydrogen_indices = [index for index, element in atom_types.items() if element == 'H']
#
## Apply constraints to non-hydrogens
#for i in range(system.getNumParticles()):
#    if i not in hydrogen_indices:
#        # Apply constraints to keep non-hydrogens fixed
#        system.addConstraint(i, i, 0.0)  # Constraint with zero distance

# Add a custom external force to fix atom positions
fix_force = openmm.CustomExternalForce('0.5*k*(x^2 + y^2 + z^2)')
fix_force.addPerParticleParameter('k')
k = 1e10 * unit.kilojoules_per_mole / unit.nanometers**2  # Large force constant
for i in range(system.getNumParticles()):
    if i not in hydrogen_indices:
#for atomIndex in non_hydrogen_atom_indices:
#      fix_force.addParticle(atomIndex, [k])
      fix_force.addParticle(i, [k])
system.addForce(fix_force)


print('Performing energy minimization...')
simulation.minimizeEnergy()

minimized_positions = simulation.context.getState(getPositions=True).getPositions()

#Now recreate the system but with no constraints
system = forcefield.createSystem(topology,
                                      nonbondedMethod=nonbondedMethod,
                                      nonbondedCutoff=nonbondedCutoff,
                                      constraints=HBonds)

print('Initiate second minimization...')
#Set the positions to the values obtained after minimization
simulation.context.setPositions(minimized_positions)

#Re-minimize but allow all atoms to move
simulation.minimizeEnergy()

#Equilibrate
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

# Simulate

print('Simulating...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
simulation.step(steps)