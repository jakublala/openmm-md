import mdtraj as md
import numpy as np
import itertools

def get_interface_contact_indices(
        pdb, 
        cutoff=0.8,
        chains='AB'
        ):
    # find the residue indices of the selected chains within a cutoff
    traj = md.load(pdb)

    assert len(chains) == 2, "Only two chains are supported"

    chain_A_indices = traj.topology.select(f'chainid {chains[0]} and name CA')
    chain_B_indices = traj.topology.select(f'chainid {chains[1]} and name CA')

    distances = md.compute_distances(traj, np.array(list(itertools.product(chain_A_indices, chain_B_indices))))
    contact_indices = np.where(distances < cutoff)

    return contact_indices
