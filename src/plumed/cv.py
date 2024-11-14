import mdtraj as md
import numpy as np
import itertools

from src.models import ContactMap, Contact, Residue

import logging
logger = logging.getLogger(__name__)

def remove_non_protein_elements(traj):
    return traj.atom_slice(traj.topology.select('protein'))

def get_contact_map(
        filename: str, 
        cutoff: float = 0.8, # in angstroms
        output_dir: str = None,
        spot1_residue_ids: list[int] = None,
        spot2_residue_ids: list[int] = None,
        ):
    """
    Creates a contact residue map for a given system.
    It outputs the indices of the contact residues.
    """
    # find the residue indices of the selected chains within a cutoff

    if output_dir is None:
        raise ValueError('Output directory is required')
    
    if spot1_residue_ids is None or spot2_residue_ids is None:
        raise ValueError('Contact indices are required')

    traj = md.load(f'{output_dir}/{filename}_solvated.pdb')
    traj = remove_non_protein_elements(traj)

    atom_pairs = list(itertools.product(spot1_residue_ids, spot2_residue_ids))
    distances = md.compute_distances(traj, atom_pairs)
    contact_atom_indices = [atom_pairs[i] for i in np.where(distances < cutoff)[1]] 

    contact_map = ContactMap()
    for i, j in contact_atom_indices:
        binder_residue_index = traj.topology.atom(i).residue.index + 1
        binder_residue_chain_id = traj.topology.atom(i).residue.chain.chain_id
        inf_residue_index = traj.topology.atom(j).residue.index + 1
        inf_residue_chain_id = traj.topology.atom(j).residue.chain.chain_id
        residue1 = Residue(
            index=binder_residue_index,
            chain_id=binder_residue_chain_id
        )
        residue2 = Residue(
            index=inf_residue_index,
            chain_id=inf_residue_chain_id
        )
        contact_map.append(
            Contact(
                residue1=residue1,
                residue2=residue2
            )
        )
    return contact_map


def get_contact_content(filename, config, output_dir):

    contact_residues = get_interface_contact_indices(
        filename, 
        cutoff=config['cutoff'], 
        output_dir=output_dir,
        spot1_indices=config['spot1_indices'],
        spot2_indices=config['spot2_indices']
        )
    
    # Generate PLUMED input content
    contact_str = '\n'.join(f"\tATOMS{i+1}=@CA-A_{a},@CA-B_{b}" for i, (a, b) in enumerate(contact_residues))

    