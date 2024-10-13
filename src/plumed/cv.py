import mdtraj as md
import numpy as np
import itertools

def get_interface_contact_indices(
        filename, 
        cutoff=0.8, # in angstroms
        chains='AB'
        ):
    # find the residue indices of the selected chains within a cutoff
    traj = md.load(f'tmp/{filename}_fixed.pdb')

    assert len(chains) == 2, "Only two chains are supported"

    chain_A_indices = traj.topology.select(f'chainid 0 and name CA')
    chain_B_indices = traj.topology.select(f'chainid 1 and name CA')

    # Ensure atom_pairs is a 2D array
    atom_pairs = np.array(list(itertools.product(chain_A_indices, chain_B_indices)))

    distances = md.compute_distances(traj, atom_pairs)
    contact_indices = atom_pairs[np.where(distances < cutoff)[1]]

    return contact_indices


