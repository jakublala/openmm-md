from pymol import cmd
import os
import numpy as np

def center_of_mass_ca(object_name, state=1):
    """
    Calculate the Center of Mass (COM) based only on the alpha carbon (CA) atoms
    of the object in a specific state of the trajectory.

    Parameters:
    - object_name: Name of the object to compute the COM for.
    - state: The state number to compute the COM for (default is 1).

    Returns:
    - COM position as a tuple (x, y, z).
    """
    # Validate the existence of the object
    if object_name not in cmd.get_names('objects'):
        raise ValueError(f"Object '{object_name}' does not exist in the current PyMOL session.")
    
    # Select only alpha carbon atoms (CA) for the specified state
    print(f"Selecting CA atoms for object: {object_name}, state: {state}")
    cmd.select('ca_atoms', f"{object_name} and name CA", state=state)
    
    # Check if any atoms were selected
    n_atoms = cmd.count_atoms('ca_atoms')
    print(f"Number of CA atoms selected: {n_atoms}")
    if n_atoms == 0:
        raise ValueError(f"No CA atoms found in object '{object_name}' for state {state}. Please check the object name and atom selection.")

    # Get the coordinates for each CA atom in the specified state
    coords = np.array(cmd.get_coords('ca_atoms', state=state))

    # Calculate the COM for CA atoms
    com = np.mean(coords, axis=0)
    print(f"Center of Mass (CA atoms) for '{object_name}', state {state}: {com}")
    return com

def get_symmetry_info(object_name):
    # Get symmetry info from _solvated.pdb object, need to first load the pdb and get it
    symmetry_info = cmd.get_symmetry(object_name)
    if symmetry_info is None:
        raise ValueError(f"Object '{object_name}' does not have symmetry information.")

    # Extract unit cell dimensions
    a, b, c = symmetry_info[0], symmetry_info[1], symmetry_info[2]

    return a, b, c

def translate_trajectory(object_name):
    """
    Translates an object across all frames in the trajectory, automatically 
    retrieving unit cell dimensions if dx, dy, or dz is not provided. It also 
    moves each frame to the COM (based on CA atoms) and then applies the translation.

    Parameters:
    - object_name: Name of the object to translate
    - dx, dy, dz: Translation values along the x, y, and z axes. If None, these
                  are computed as half the unit cell dimensions.
    """
    a, b, c = get_symmetry_info(f"{object_name}_solvated")

    dx = a / 2
    dy = b / 2
    dz = c / 2

    print(f"Translating '{object_name}' by: dx={dx}, dy={dy}, dz={dz}")

    # Apply translation to all frames
    num_states = cmd.count_states(object_name)
    for state in range(1, num_states + 1):  # PyMOL states are 1-indexed
        # Calculate the center of mass for the current state
        dx_, dy_, dz_ = center_of_mass_ca(object_name, state)
        
        # Subtract center of mass from current state
        cmd.alter_state(state, object_name,
                        f"(x, y, z) = (x - {dx_}, y - {dy_}, z - {dz_})")
        
        # Move to center of unit cell
        cmd.alter_state(state, object_name,
                        f"(x, y, z) = (x + {dx}, y + {dy}, z + {dz})")
    cmd.refresh()
    print(f"Translation complete for {num_states} frames of '{object_name}'.")

def replicate_trajectory(object_name):
    """
    Replicates the given object in a 3x3x3 grid using its periodic boundary conditions
    across all frames in the trajectory.

    Parameters:
    - object_name: Name of the object to replicate.
    """
    a, b, c = get_symmetry_info(f"{object_name}_solvated")
    
    # Create a 3x3x3 grid of the object
    for i in range(-1, 2):  # -1, 0, 1
        for j in range(-1, 2):  # -1, 0, 1
            for k in range(-1, 2):  # -1, 0, 1
                if i == 0 and j == 0 and k == 0:
                    continue  # Skip the original object position
                new_object_name = f"{object_name}_copy_{i+1}_{j+1}_{k+1}"
                
                cmd.create(new_object_name, object_name, 0, 0)
                
                # Changed: Apply translation to all states at once using state=0
                # Removed: loop over individual states
                translation = f"(x + {i * a}, y + {j * b}, z + {k * c})"
                cmd.alter_state(0, new_object_name, f"(x, y, z) = {translation}")

    print(f"Replicated '{object_name}' into a 3x3x3 grid across all frames.")

import re

def extract_residues(file_content):
    # Extract binding site residues
    binding_pattern = r"select binding_site, resi ([\d+]+) and chain"
    binding_match = re.search(binding_pattern, file_content)
    binder_spots = binding_match.group(1).split('+') if binding_match else []
    
    # Extract binder contact residues
    contact_pattern = r"select binder_contacts, resi ([\d+]+) and chain"
    contact_match = re.search(contact_pattern, file_content)
    target_spots = contact_match.group(1).split('+') if contact_match else []
    
    return binder_spots, target_spots

def colour_object(object_name, commands_file):
    # open *_commands.txt file
    with open(commands_file, 'r') as file:
        commands = file.read()
        binder_spots, target_spots = extract_residues(commands)

    # Join residue numbers with + to create a valid PyMOL selection
    binder_spots_str = '+'.join(binder_spots)
    target_spots_str = '+'.join(target_spots)

    # compute the number of different chains
    num_chains = len(set(cmd.get_chains(object_name)))

    if num_chains == 2:
        cmd.color("blue", f"{object_name} and chain A")
        cmd.color("green", f"{object_name} and chain B")
        cmd.color("magenta", f"{object_name} and resi {binder_spots_str} and chain A")
        cmd.color("yellow", f"{object_name} and resi {target_spots_str} and chain B")
    elif num_chains == 1:
        cmd.color("blue", f"{object_name} and chain A")
        cmd.color("magenta", f"{object_name} and resi {binder_spots_str} and chain A")
        cmd.color("yellow", f"{object_name} and resi {target_spots_str} and chain A")
    else:
        raise ValueError(f"Expected 1 or 2 chains, got {num_chains}")
    

def load_trajectory(replicate=False):
    """
    Automatically loads the PDB and DCD files from the current directory, centers
    the protein based on the CA atoms, translates the trajectory to the center of
    the unit cell, and executes necessary PyMOL commands.
    """
    # Get the current directory and look for PDB and DCD files
    current_directory = os.getcwd()
    pdb_file = None
    dcd_file = None

    # Search for PDB and DCD files in the current directory
    for file in os.listdir(current_directory):
        if file.endswith("fixed.pdb"):
            pdb_file = file
        elif file.endswith("solvated.pdb"):
            solvated_pdb_file = file
        elif file.endswith(".dcd"):
            dcd_file = file
        elif file.endswith("commands.txt"):
            commands_file = file

    # Ensure both PDB and DCD files are found
    if pdb_file is None or dcd_file is None:
        raise FileNotFoundError("Both a PDB and a DCD file must be present in the current directory.")

    print(f"Loading PDB: {pdb_file} and DCD: {dcd_file}")

    abs_path = os.path.abspath(current_directory)
    path_parts = abs_path.split(os.sep)
    system_name = path_parts[-2]
    # Load the PDB and DCD files
    cmd.load(pdb_file, system_name)
    cmd.load_traj(dcd_file, system_name)

    # Load the solvated pdb file and hide it
    cmd.load(solvated_pdb_file, f"{system_name}_solvated")
    cmd.hide("everything", f"{system_name}_solvated")

    # Reset the view and show the unit cell
    cmd.reset()
    cmd.show("cell", f"{system_name}_solvated")

    colour_object(system_name, commands_file)

    # Center the trajectory to the COM based on CA atoms
    translate_trajectory(system_name)

    print("Trajectory loaded, centered, and translated.")

    if replicate:
        replicate_trajectory(system_name)