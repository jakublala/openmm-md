import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix
import numpy as np
import logging
from typing import Optional

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

def get_contact_map(
        filename: str, 
        cutoff: float = 0.8, # in angstroms
        output_dir: str = None,
        spot1_residues: Optional[Segment] = None,
        spot2_residues: Optional[Segment] = None,
        ):
    """
    Creates a contact residue map for a given system.
    It outputs the indices of the contact residues.
    """
    # find the residue indices of the selected chains within a cutoff

    if output_dir is None:
        raise ValueError('Output directory is required')

    universe = mda.Universe(f'{output_dir}/{filename}_solvated.pdb')
    chain_ids = np.concatenate(universe.segments.chainIDs)
    universe = get_CA_universe(universe)

    cmatrix = contact_matrix(
        universe.atoms,
        cutoff=cutoff*10,
        box=universe.dimensions, 
        )
    np.fill_diagonal(cmatrix, False)
    contact_residue_ids = np.argwhere(cmatrix)
    contact_residue_ids += 1 # mda is 1-indexed
    if spot1_residues is None and spot2_residues is None:
        logger.info('No spot definition provided, using all contacts')
    else:
        contact_residue_ids = filter_contact_residue_ids(contact_residue_ids, spot1_residues, spot2_residues)

    contact_map = ContactMap()
    for spot1_residue_index, spot2_residue_index in contact_residue_ids:
        spot1_residue_chain_id = chain_ids[spot1_residue_index]
        spot2_residue_chain_id = chain_ids[spot2_residue_index]
        assert spot1_residue_chain_id is not None, f"Chain ID is None for spot1 residue {spot1_residue_index}"
        assert spot2_residue_chain_id is not None, f"Chain ID is None for spot2 residue {spot2_residue_index}"
        residue1 = Residue(
            index=spot1_residue_index,
            chain_id=spot1_residue_chain_id,
            indexing=1
        )
        residue2 = Residue(
            index=spot2_residue_index,
            chain_id=spot2_residue_chain_id,
            indexing=1
        )
        contact_map.contacts.append(
            Contact(
                residue1=residue1,
                residue2=residue2
            )
        )
    return contact_map
