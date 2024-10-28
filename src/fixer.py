from pdbfixer import PDBFixer
from openmm.app import PDBFile
import fire
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import os

def fixer(filepath=None, output_dir=None):
    if filepath is None:
        raise ValueError('filepath is required')
    
    file = os.path.basename(filepath).split('.')[0]

    print('Splitting chains...')
    structure = PDBParser().get_structure('protein', filepath)
    
    all_ids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                all_ids.append(residue.id[1])
    
    indices_to_break = []
    for i in range(len(all_ids) - 1):
        if abs(all_ids[i] - all_ids[i + 1]) > 2:
            indices_to_break.append(i)

    chains = {'A': [], 'B': [], 'C': []}

    new_structure = Structure(file)
    new_model = Model(0)
    new_structure.add(new_model)

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

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(f"{output_dir}/{file}_fixed.pdb")
    fixer = PDBFixer(filename=f"{output_dir}/{file}_fixed.pdb")

    print('Fixing C-terminus...')
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in keys:
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]
    
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()

    with open(f'{output_dir}/{file}_fixed.pdb', 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


if __name__ == '__main__':
    fire.Fire(fixer)