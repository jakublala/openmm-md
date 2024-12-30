import openmm
from openmm import app, unit
import copy
import math
import logging
import numpy as np
from openmmtools.constants import kB
import openmmtools
from openmmtools import cache, mcmc, multistate
from openmmtools.multistate import ReplicaExchangeSampler
from openmmtools.states import GlobalParameterState, SamplerState, ThermodynamicState, CompoundThermodynamicState
from openmm.app import StateDataReporter, PDBReporter
import sys
import os
import time

"""
This is implemented from the REST tutorial on the OpenMM website.
https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/Running_a_REST_simulation.html
"""

class RESTState(GlobalParameterState):
    lambda_rest_bonds = GlobalParameterState.GlobalParameter('lambda_rest_bonds', standard_value=1.0)
    lambda_rest_angles = GlobalParameterState.GlobalParameter('lambda_rest_angles', standard_value=1.0)
    lambda_rest_torsions = GlobalParameterState.GlobalParameter('lambda_rest_torsions', standard_value=1.0)
    lambda_rest_electrostatics = GlobalParameterState.GlobalParameter('lambda_rest_electrostatics', standard_value=0.0)
    lambda_rest_sterics = GlobalParameterState.GlobalParameter('lambda_rest_sterics', standard_value=0.0)

    def set_rest_parameters(self, beta_m, beta_0):
        """Set all defined lambda parameters to the given value.

        The undefined parameters (i.e. those being set to None) remain undefined.

        Parameters
        ----------
        new_value : float
            The new value for all defined parameters.
        """
        lambda_functions = {'lambda_rest_bonds': lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0),
                 'lambda_rest_angles' : lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0),
                 'lambda_rest_torsions' : lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0),
                 'lambda_rest_electrostatics' : lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0) - 1,
                 'lambda_rest_sterics' : lambda beta_m, beta_0 : beta_m / beta_0 - 1
                 }

        for parameter_name in self._parameters:
            if self._parameters[parameter_name] is not None:
                new_value = lambda_functions[parameter_name](beta_m, beta_0)
                setattr(self, parameter_name, new_value)

def get_rest_identifier(atoms, rest_atoms):
    """
    For a given atom or set of atoms, get the rest_id which is a list of binary ints that defines which
    (mutually exclusive) set the atom(s) belong to.

    If there is a single atom, the sets are: [is_rest, is_nonrest]
    If there is a set of atoms, the sets are: [is_rest, is_inter, is_nonrest]
    
    is_rest means we are in the solute.
    is_nonrest means we are in the solvent.
    is_inter means we are in the interface, i.e. the solute-solvent interaction, that's why we need a set of atoms.

    Example: if there is a single atom that is in the nonrest set, the rest_id is [0, 1]

    Arguments
    ---------
    atoms : set or int
        a set of hybrid atom indices or single atom
    rest_atoms : set or list
        a list (or list-like) of atoms whose interactions will be scaled by REST
    Returns
    -------
    rest_id : list
        list of binaries indicating which set the atom(s) belong to
    """

    if isinstance(atoms, int):
        rest_id = [0, 1] # Set the default rest_id to non-REST
        if atoms in rest_atoms:
            rest_id = [1, 0]
        return rest_id

    elif isinstance(atoms, set):
        rest_id = [0, 0, 1] # Set the default rest_id to non-REST
        if atoms.intersection(rest_atoms) != set(): # At least one of the atoms is REST
            if atoms.issubset(rest_atoms): # All atoms are REST
                rest_id = [1, 0, 0]
            else: # At least one (but not all) of the atoms is are REST
                rest_id = [0, 1, 0]
        return rest_id

    else:
        raise Exception(f"atoms is of type {type(atoms)}, but only `int` and `set` are allowable")

import openmm as mm

def apply_scaling(
        system: mm.System, 
        rest_atoms: list, # list of solute atom indices
        ):
    # barostat = openmm.MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 50)
    # system.addForce(barostat)

    print(system.getForces())

    print('len(rest_atoms)', len(rest_atoms))

    # Create REST system
    rest_system = openmm.System()

    # Create dict of vanilla system forces (for easy retrieval of force objects later)
    system_forces = {type(force).__name__ : force for force in system.getForces()}

    # Add particles
    for particle_idx in range(system.getNumParticles()):
        particle_mass = system.getParticleMass(particle_idx)
        rest_system.addParticle(particle_mass)

    # Copy barostat
    if "MonteCarloBarostat" in system_forces:
        barostat = copy.deepcopy(system_forces["MonteCarloBarostat"])
        rest_system.addForce(barostat)

    # No need to add thermostat, as it's the integrator?

    # Copy box vectors
    box_vectors = system.getDefaultPeriodicBoxVectors()
    rest_system.setDefaultPeriodicBoxVectors(*box_vectors)

    # Copy constraints
    for constraint_idx in range(system.getNumConstraints()):
        atom1, atom2, length = system.getConstraintParameters(constraint_idx)
        rest_system.addConstraint(atom1, atom2, length)
        
    # Copy PLUMED force
    if "PlumedForce" in system_forces:
        plumed_force = copy.deepcopy(system_forces["PlumedForce"])
        rest_system.addForce(plumed_force)

    # Define the custom expression
    bond_expression = "rest_scale * (K / 2) * (r - length)^2;"
    bond_expression += "rest_scale = is_rest * lambda_rest_bonds * lambda_rest_bonds " \
                    "+ is_inter * lambda_rest_bonds " \
                    "+ is_nonrest;"

    # Create custom force
    rest_bond_force = openmm.CustomBondForce(bond_expression)
    rest_system.addForce(rest_bond_force)

    # Add global parameters
    rest_bond_force.addGlobalParameter("lambda_rest_bonds", 1.0)

    # Add per-bond parameters for rest scaling
    rest_bond_force.addPerBondParameter("is_rest")
    rest_bond_force.addPerBondParameter("is_inter")
    rest_bond_force.addPerBondParameter("is_nonrest")

    # Add per-bond parameters for defining bond energy
    rest_bond_force.addPerBondParameter('length')  # equilibrium bond length
    rest_bond_force.addPerBondParameter('K')  # force constant

    # Get vanilla system bond force
    bond_force = system_forces['HarmonicBondForce']

    # Set periodicity
    if bond_force.usesPeriodicBoundaryConditions():
        rest_bond_force.setUsesPeriodicBoundaryConditions(True)

    # Add bonds to rest_system
    for term_idx in range(bond_force.getNumBonds()):
        # Get the bond parameters and rest id
        p1, p2, r0, k = bond_force.getBondParameters(term_idx)
        idx_set = set([p1, p2])
        rest_id = get_rest_identifier(idx_set, rest_atoms)

        # Add the bond
        bond_term = (p1, p2, rest_id + [r0, k])
        rest_bond_force.addBond(*bond_term)
        
    # Define the custom expression
    angle_expression = "rest_scale * (K / 2) * (theta - theta0)^2;"
    angle_expression += "rest_scale = is_rest * lambda_rest_angles * lambda_rest_angles " \
                        "+ is_inter * lambda_rest_angles " \
                        "+ is_nonrest;"

    # Create custom force
    rest_angle_force = openmm.CustomAngleForce(angle_expression)
    rest_system.addForce(rest_angle_force)

    # Add global parameters
    rest_angle_force.addGlobalParameter("lambda_rest_angles", 1.0)

    # Add per-angle parameters for rest scaling
    rest_angle_force.addPerAngleParameter("is_rest")
    rest_angle_force.addPerAngleParameter("is_inter")
    rest_angle_force.addPerAngleParameter("is_nonrest")

    # Add per-angle parameters for defining angle energy
    rest_angle_force.addPerAngleParameter('theta0')  # equilibrium angle
    rest_angle_force.addPerAngleParameter('K')  # force constant

    # Get vanilla system angle force
    angle_force = system_forces['HarmonicAngleForce']

    # Set periodicity
    if angle_force.usesPeriodicBoundaryConditions():
        rest_angle_force.setUsesPeriodicBoundaryConditions(True)

    # Add angles to rest_system
    for term_idx in range(angle_force.getNumAngles()):
        # Get the angle parameters and rest id
        p1, p2, p3, theta0, k = angle_force.getAngleParameters(term_idx)
        idx_set = set([p1, p2, p3])
        rest_id = get_rest_identifier(idx_set, rest_atoms)

        # Add the angle
        angle_term = (p1, p2, p3, rest_id + [theta0, k])
        rest_angle_force.addAngle(*angle_term)
        
    # Define the custom expression
    torsion_expression = "rest_scale * U;"
    torsion_expression += "rest_scale = is_rest * lambda_rest_torsions * lambda_rest_torsions " \
                        "+ is_inter * lambda_rest_torsions " \
                        "+ is_nonrest;"
    torsion_expression += "U = (K * (1 + cos(periodicity * theta - phase)));"

    # Create custom force
    rest_torsion_force = openmm.CustomTorsionForce(torsion_expression)
    rest_system.addForce(rest_torsion_force)

    # Add global parameters
    rest_torsion_force.addGlobalParameter("lambda_rest_torsions", 1.0)

    # Add per-torsion parameters for rest scaling
    rest_torsion_force.addPerTorsionParameter("is_rest")
    rest_torsion_force.addPerTorsionParameter("is_inter")
    rest_torsion_force.addPerTorsionParameter("is_nonrest")

    # Add per-torsion parameters for defining torsion energy
    rest_torsion_force.addPerTorsionParameter('periodicity')
    rest_torsion_force.addPerTorsionParameter('phase') # phase offset
    rest_torsion_force.addPerTorsionParameter('K') # force constant

    # Get vanilla system torsion force
    torsion_force = system_forces['PeriodicTorsionForce']

    # Set periodicity
    if torsion_force.usesPeriodicBoundaryConditions():
        rest_torsion_force.setUsesPeriodicBoundaryConditions(True)

    # Add torsions to rest_system
    for torsion_idx in range(torsion_force.getNumTorsions()):
        # Get the torsion parameters and rest id
        p1, p2, p3, p4, periodicity, phase, K = torsion_force.getTorsionParameters(torsion_idx)
        idx_set = set([p1, p2, p3, p4])
        rest_id = get_rest_identifier(idx_set, rest_atoms)

        # Add torsion
        torsion_term = (p1, p2, p3, p4, rest_id + [periodicity, phase, K])
        rest_torsion_force.addTorsion(*torsion_term)
        
    # Create nonbonded force
    rest_nonbonded_force = openmm.NonbondedForce()
    rest_system.addForce(rest_nonbonded_force)

    # Get vanilla system nonbonded force
    nonbonded_force = system_forces['NonbondedForce']

    # Set the nonbonded method and related parameters
    nonbonded_method = nonbonded_force.getNonbondedMethod()
    rest_nonbonded_force.setNonbondedMethod(nonbonded_method)
    if nonbonded_method != openmm.NonbondedForce.NoCutoff:
        epsilon_solvent = nonbonded_force.getReactionFieldDielectric()
        cutoff = nonbonded_force.getCutoffDistance()
        rest_nonbonded_force.setReactionFieldDielectric(epsilon_solvent)
        rest_nonbonded_force.setCutoffDistance(cutoff)
    if nonbonded_method in [openmm.NonbondedForce.PME, openmm.NonbondedForce.Ewald]:
        [alpha_ewald, nx, ny, nz] = nonbonded_force.getPMEParameters()
        delta = nonbonded_force.getEwaldErrorTolerance()
        rest_nonbonded_force.setPMEParameters(alpha_ewald, nx, ny, nz)
        rest_nonbonded_force.setEwaldErrorTolerance(delta)

    # Copy switching function from vanilla system
    switch_bool = nonbonded_force.getUseSwitchingFunction()
    rest_nonbonded_force.setUseSwitchingFunction(switch_bool)
    if switch_bool:
        switching_distance = nonbonded_force.getSwitchingDistance()
        rest_nonbonded_force.setSwitchingDistance(switching_distance)

    # Copy dispersion correction
    dispersion_bool = nonbonded_force.getUseDispersionCorrection()
    rest_nonbonded_force.setUseDispersionCorrection(dispersion_bool)

    # Add global parameters
    rest_nonbonded_force.addGlobalParameter('lambda_rest_electrostatics', 0.)
    rest_nonbonded_force.addGlobalParameter('lambda_rest_sterics', 0.)


    # Add nonbondeds to rest_system
    for particle_idx in range(nonbonded_force.getNumParticles()):
        # Get the nonbonded parameters and rest id
        q, sigma, epsilon = nonbonded_force.getParticleParameters(particle_idx)
        rest_id = get_rest_identifier(particle_idx, rest_atoms)

        # Add particles and offsets
        if rest_id == [0, 1]: # nonrest
            rest_nonbonded_force.addParticle(q, sigma, epsilon)

        else: # rest
            rest_nonbonded_force.addParticle(q, sigma, epsilon)
            rest_nonbonded_force.addParticleParameterOffset('lambda_rest_electrostatics', particle_idx, q, 0.0*sigma, epsilon*0.0)
            rest_nonbonded_force.addParticleParameterOffset('lambda_rest_sterics', particle_idx, q*0.0, 0.0*sigma, epsilon)

    # Handle exceptions
    for exception_idx in range(nonbonded_force.getNumExceptions()):
        # Get exception parameters and rest id
        p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(exception_idx)
        idx_set = set([p1, p2])
        rest_id = get_rest_identifier(idx_set, rest_atoms)

        # Add exceptions and offsets
        exc_idx = rest_nonbonded_force.addException(p1, p2, chargeProd, sigma, epsilon)
        if rest_id == [0, 0, 1]: # nonrest
            pass

        elif rest_id == [1, 0, 0]: # rest
            rest_nonbonded_force.addExceptionParameterOffset('lambda_rest_sterics', exc_idx, chargeProd, 0.0*sigma, epsilon)

        elif rest_id == [0, 1, 0]: # inter
            rest_nonbonded_force.addExceptionParameterOffset('lambda_rest_electrostatics', exc_idx, chargeProd, 0.0*sigma, epsilon)
    

    return rest_system


import tempfile
from openmm.app import PDBxFile
import MDAnalysis as mda
import re
def parse_plumed_cmap(plumed_input_file):
    with open(plumed_input_file, 'r') as file:
        script = file.read()

    chainA_residue_ids = set(map(int, re.findall(r'@CA-A_(\d+)', script)))
    chainB_residue_ids = set(map(int, re.findall(r'@CA-B_(\d+)', script)))
    
    # these are local residue indices, not global residue indices
    return chainA_residue_ids, chainB_residue_ids
    

def get_rest_atoms(system, pdb, plumed_input_file):
    """
        Heat up the atoms that are part of the CMAP.
        Assume chain A is the binder; then chain B is the target.
        We heat up entire chain A, but for chain B only the atoms that are part of the CMAP.
        We do this awkwardly by parsing the PLUMED input file.
    """
    # HACK: hardcoded path to the PLUMED input file
    _, chainB_residue_ids = parse_plumed_cmap(plumed_input_file)
    # ignore chain A as we heat up the entire binder
    universe = mda.Universe(pdb)
    # get all chain A atoms
    chainA_atoms = universe.select_atoms('chainid A').indices.tolist()
    # get all chain B atoms that are part of the CMAP
    chainB_atoms = universe.select_atoms(f'chainid B and resid {" ".join(str(i) for i in chainB_residue_ids)}').indices.tolist()
    rest_atoms = chainA_atoms + chainB_atoms
    # these are 0-indexed atom indices
    return rest_atoms
    
from openmmplumed import PlumedForce
def get_thermodynamic_and_sampler_states(system, temperatures, equilibriated_file, plumed_input_file):
    
    # HACK: assuming that the temperatures are in ascending order
    T_min = temperatures[0]

    beta_0 = 1/(kB*T_min)
    thermodynamic_state_list = []
    sampler_state_list = []

    pdb = PDBxFile(equilibriated_file)

    rest_atoms = get_rest_atoms(system, pdb, plumed_input_file)
    system = apply_scaling(system, rest_atoms)

    rest_state = RESTState.from_system(system)
    thermostate = ThermodynamicState(system, temperature=T_min)
    compound_thermodynamic_state = CompoundThermodynamicState(thermostate, composable_states=[rest_state])
    sampler_state =  SamplerState(pdb.positions, box_vectors=system.getDefaultPeriodicBoxVectors())

    plumed_input_file = plumed_input_file.replace('.dat', '')
    for i, temperature in enumerate(temperatures):
        # Create a thermodynamic state with REST interactions scaled to the given temperature
        beta_m = 1/(kB*temperature)
        compound_thermodynamic_state_copy = copy.deepcopy(compound_thermodynamic_state)
        compound_thermodynamic_state_copy.set_rest_parameters(beta_m, beta_0)

        thermodynamic_state_list.append(compound_thermodynamic_state_copy)

        # Generate a sampler_state with minimized positions
        context, integrator = cache.global_context_cache.get_context(compound_thermodynamic_state_copy)
        sampler_state.apply_to_context(context, ignore_velocities=True)
        openmm.LocalEnergyMinimizer.minimize(context)
        sampler_state.update_from_context(context)
        sampler_state_list.append(copy.deepcopy(sampler_state))

        # add PLUMED force to the last thermodynamic state's system
        system = thermodynamic_state_list[-1].get_system()
        with open(f'{plumed_input_file}_{i}.dat', 'r') as file:
            plumed_script = file.read()
        system.addForce(PlumedForce(plumed_script))
        thermodynamic_state_list[-1].set_system(system)

    return thermodynamic_state_list, sampler_state_list

