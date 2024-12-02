from pdbfixer import PDBFixer
from openmm.app import PDBFile
import fire
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import os

import logging
logger = logging.getLogger(__name__)

# Add these lines to configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def get_concat_sequence_from_structure(structure):
    # biopython
    sequence = ''
    for model in structure:
        for chain in model:
            for residue in chain:
                sequence += residue.get_resname()
    return sequence

def get_concat_sequence_from_topology(topology):
    # openmm
    sequence = ''
    for chain in topology.chains():
        for residue in chain.residues():
            sequence += residue.name
    return sequence

def assert_sequence_identity(original_sequence, fixed_sequence):
    if original_sequence != fixed_sequence:
        raise ValueError('Fixer.py: Fixing the PDB introduced / deleted some residues!')

def get_all_ids(structure):
    all_ids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                all_ids.append(residue.id[1])
    return all_ids

def get_break_indices(structure):
    all_ids = get_all_ids(structure)
    indices_to_break = []
    for i in range(len(all_ids) - 1):
        if abs(all_ids[i] - all_ids[i + 1]) > 2:
            indices_to_break.append(i)
    return indices_to_break



def check_already_fixed_and_split(structure):
    break_indices = get_break_indices(structure)
    if len(break_indices) == 0:
        logger.info('No chain break indices found, assuming chains are already split and correct.')
        return True
    elif len(break_indices) == 1:
        logger.info('One chain break index found, assuming chains are split and correct.')
        return True
    else:
        raise ValueError('Fixer.py: More than one chain break indices found!')


def fixer(filepath=None, output_dir=None, split_chains=True):
    if filepath is None:
        raise ValueError('filepath is required')
    if output_dir is None:
        raise ValueError('output_dir is required')
    
    file = os.path.basename(filepath).split('.')[0]

    logger.info(f"Fixing PDB protein file {file}...")

    structure = PDBParser().get_structure('protein', filepath)

    original_sequence = get_concat_sequence_from_structure(structure)

    # # before splitting chains, check whether the chains are actually not already split and correct
    # if check_already_fixed_and_split(structure):
    #     logger.info('Chains are already split and correct.')
    #     split_chains = False

    if split_chains:
        logger.info('Splitting chains... (usually run from ESMFold outputs where chains are not split)')

        indices_to_break = get_break_indices(structure)
        chains = {'A': [], 'B': [], 'C': []}

        new_structure = Structure(file)
        new_model = Model(0)
        new_structure.add(new_model)

        for model in structure:
            for index, chain in enumerate(model):
                print(chain)
                if index != 0:
                    raise ValueError('Fixer.py: More than one model in the PDB!')

        for model in structure:
            for chain in model:
                residues = [i for i in chain.get_residues()]
                if len(indices_to_break) > 1:
                    chains['A'] = residues[0:indices_to_break[0]+1]
                    chains['B'] = residues[indices_to_break[0]+1:indices_to_break[1]+1]
                    chains['C'] = residues[indices_to_break[1]+1:]
                else:
                    chains['A'] = residues[0:indices_to_break[0]+1]
                    chains['B'] = residues[indices_to_break[0]+1:]

        residue_id = 1
        for chain_id, residues in chains.items():
            new_chain = Chain(chain_id)
            for residue in residues:
                residue.id = (residue.id[0], residue_id, residue.id[2])
                new_chain.add(residue)
                residue_id += 1
            new_model.add(new_chain)
    else:
        logger.info('Not splitting chains, assuming that this is a single construct of a single chain; or everything is already split correctly.')
        new_structure = structure

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(f"{output_dir}/{file}_fixed.pdb")
    fixer = PDBFixer(filename=f"{output_dir}/{file}_fixed.pdb")
    logger.info('Fixing C-terminus...')
    logger.info('Finding missing residues...')
    fixer.findMissingResidues()
    
    logger.info('Getting chains and keys...')
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    
    logger.info('Processing missing residues...')
    for key in keys:
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]
    
    logger.info('Finding missing atoms...')
    fixer.findMissingAtoms()
    
    logger.info('Adding missing atoms...')
    fixer.addMissingAtoms()
    
    logger.info('Adding missing hydrogens...')
    fixer.addMissingHydrogens()

    fixed_sequence = get_concat_sequence_from_topology(fixer.topology)
    assert_sequence_identity(original_sequence, fixed_sequence)



    logger.info('Writing fixed PDB file...')
    with open(f'{output_dir}/{file}_fixed.pdb', 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    logger.info("Fixed PDB successfully written.")


if __name__ == '__main__':
    fire.Fire(fixer)