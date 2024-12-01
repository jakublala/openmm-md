import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix
import numpy as np
import logging
from typing import Optional, Literal

import os
from src.models import ContactMap, Contact, Residue
from src.models import Segment

logger = logging.getLogger(__name__)

def get_CA_universe(universe):
    # note: there's a bug in mdtraj where traj.atom_slice(traj.topology.select('protein'))
    # loses the information about the chain id
    # hence moved to MDAnalysis
    return mda.Merge(universe.select_atoms('name CA'))

def filter_contact_residue_ids(
        contact_residue_ids: np.ndarray, 
        spot1_residues: Segment, 
        spot2_residues: Segment
        ) -> np.ndarray:
    mask1 = np.isin(contact_residue_ids[:, 0], spot1_residues.ids)
    mask2 = np.isin(contact_residue_ids[:, 1], spot2_residues.ids)
    return contact_residue_ids[mask1 & mask2]

def filter_interchain_contacts(
        contact_residue_ids: np.ndarray, 
        chain_ids: np.ndarray
        ) -> np.ndarray:
    contact_residue_ids_return = []
    for c1, c2 in contact_residue_ids:
        # turn to 0-indexed
        if chain_ids[c1 - 1] != chain_ids[c2 - 1]:
            contact_residue_ids_return.append((c1, c2))
    return np.array(contact_residue_ids_return)
    

def get_contact_map(
        filename: str, 
        cutoff: float, # in nm
        output_dir: str = None,
        spot1_residues: Optional[Segment] = None,
        spot2_residues: Optional[Segment] = None,
        mode: Literal['single-chain', 'two-chain'] = None,
        ):
    """
    Creates a contact residue map for a given system.
    It outputs the indices of the contact residues.
    """
    # find the residue indices of the selected chains within a cutoff

    if output_dir is None:
        raise ValueError('Output directory is required')
    
    if mode is None:
        raise ValueError('Mode is required')
    
    if os.path.exists(f'{output_dir}/{filename}_equilibrated.pdb'):
        logger.info('Using equilibrated pdb file to compute contact map')
        universe = mda.Universe(f'{output_dir}/{filename}_equilibrated.pdb')
    else:
        raise ValueError(f'Cannot compute contact map with no equilibrated pdb file found for {filename}')

    full_universe = universe
    CA_universe = get_CA_universe(universe)
    global_to_local_map = CA_universe.residues.resids
    chain_ids = np.concatenate(CA_universe.segments.chainIDs)
    
    include_cutoff = 1.0*cutoff*10 # TODO: maybe include more contacts
    cmatrix = contact_matrix(
        CA_universe.atoms,
        cutoff=include_cutoff,
        box=CA_universe.dimensions, 
        )
    np.fill_diagonal(cmatrix, False)
    cmatrix = np.triu(cmatrix)
    contact_residue_ids = np.argwhere(cmatrix)
    contact_residue_ids += 1 # mda is 1-indexed

    for c1, c2 in contact_residue_ids:
        assert c1 < c2, f"Contact residue indices are not in the correct order: {c1} < {c2}"
    if spot1_residues is None and spot2_residues is None:
        logger.info('No spot definition provided, using all contacts')
    else:
        contact_residue_ids = filter_contact_residue_ids(contact_residue_ids, spot1_residues, spot2_residues)
    if mode == 'two-chain':
        contact_residue_ids = filter_interchain_contacts(contact_residue_ids, chain_ids)

    contact_map = ContactMap()
    for spot1_residue_index, spot2_residue_index in contact_residue_ids:
        spot1_residue_chain_id = chain_ids[spot1_residue_index - 1]
        spot2_residue_chain_id = chain_ids[spot2_residue_index - 1]
        assert spot1_residue_chain_id is not None, f"Chain ID is None for spot1 residue {spot1_residue_index}"
        assert spot2_residue_chain_id is not None, f"Chain ID is None for spot2 residue {spot2_residue_index}"
        local_index1 = global_to_local_map[spot1_residue_index - 1]
        local_index2 = global_to_local_map[spot2_residue_index - 1]
        residue1 = Residue(
            index=local_index1,
            global_index=spot1_residue_index,
            chain_id=spot1_residue_chain_id,
            indexing=1,
            atom_indices=full_universe.select_atoms(f'chainid {spot1_residue_chain_id} and resid {local_index1}').indices
        )
        residue2 = Residue(
            index=local_index2,
            global_index=spot2_residue_index,
            chain_id=spot2_residue_chain_id,
            indexing=1,
            atom_indices=full_universe.select_atoms(f'chainid {spot2_residue_chain_id} and resid {local_index2}').indices
        )
        contact_map.contacts.append(
            Contact(
                residue1=residue1,
                residue2=residue2
            )
        )
    return contact_map