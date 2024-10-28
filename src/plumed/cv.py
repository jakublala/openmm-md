import mdtraj as md
import numpy as np
import itertools
import os

def get_interface_contact_indices(
        filename, 
        cutoff=0.8, # in angstroms
        chains='AB',
        output_dir=None
        ):
    # find the residue indices of the selected chains within a cutoff

    if output_dir is None:
        raise ValueError('Output directory is required')

    # HACK: as of now, mdtraj cannot deal with overflow residue numbers
    # so for this task, let's manually remove all non-protein elements
    import tempfile

    # go through all lines of the pdb until you find TER of chain B and then remove all lines until you find the next chain A TER and add END
    with open(f'{output_dir}/{filename}_solvated.pdb', 'r') as file:
        lines = file.readlines()

    # find the index of the TER line for chain B
    ter_b_index = None
    for i, line in enumerate(lines):
        if line.startswith('TER') and line.split()[3] == chains[1]:
            ter_b_index = i
            break
    lines = lines[:ter_b_index]
    
    with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as temp_file:
        for line in lines:
            temp_file.write(line.encode())
        temp_file.write(b'END\n')
        temp_filename = temp_file.name
    
    traj = md.load(temp_filename)

    assert len(chains) == 2, "Only two chains are supported"

    chain_A_indices = traj.topology.select(f'chainid 0 and name CA')
    chain_B_indices = traj.topology.select(f'chainid 1 and name CA')

    # Ensure atom_pairs is a 2D array
    atom_pairs = np.array(list(itertools.product(chain_A_indices, chain_B_indices)))

    distances = md.compute_distances(traj, atom_pairs)
    contact_indices = atom_pairs[np.where(distances < cutoff)[1]]

    os.remove(temp_filename)

    return contact_indices


