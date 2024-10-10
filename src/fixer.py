from pdbfixer import PDBFixer
from openmm.app import PDBFile
import fire
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


def fixer(
        filename=None
):
    if filename is None:
        raise ValueError('Filename is required')


    # only run if AB in filename, since that means we have two chains
    if 'AB' in filename:
        print('Splitting chains...')
        # SEE https://stackoverflow.com/questions/74735845/splitting-and-renaming-protein-chain-with-biopythons-biopdb
        structure = PDBParser().get_structure(filename, f'input/{filename}.pdb')
        # get residue id where to break
        last_id = 0
        res_to_chain_B = None
        chainB = False
        for model in structure:
            for chain in model:
                for residue in chain:
                    if chainB:
                        res_to_chain_B.append(residue)
                        continue
                    elif abs(residue.id[1] - last_id) > 2:
                        chainB = True
                        res_to_chain_B = [residue]
                    last_id = residue.id[1]

        ### SEE https://stackoverflow.com/questions/25884758/deleteing-residue-from-pdb-using-biopython-library
        for model in structure:
            for chain in model:
                [chain.detach_child(res.get_id()) for res in res_to_chain_B]

        ### SEE https://stackoverflow.com/questions/33364370/how-to-add-chain-id-in-pdb
        my_chain = Chain("B")
        model.add(my_chain)
        for index, res in enumerate(res_to_chain_B):
            # reindex
            res.id = (' ', index+1, ' ')
            my_chain.add(res)
        io = PDBIO()
        io.set_structure(model)
        io.save(f"tmp/{filename}_fixed.pdb")
        fixer = PDBFixer(filename=f"tmp/{filename}_fixed.pdb")
    else:
        fixer = PDBFixer(filename=f'input/{filename}.pdb')

    # this script added the C-terminus oxygen atom
    print('Fixing C-terminus...')
    fixer.findMissingResidues()
    # only deal with missing residues in the middle of the chain, not at the start or end
    # if this is removed, multiple chains become merged into one chain
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in keys:
        print(key)
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            print(f"Removing missing residue at chain {key[0]} and residue {key[1]}")
            del fixer.missingResidues[key]
    
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()

    with open(f'tmp/{filename}_fixed.pdb', 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


if __name__ == '__main__':
    fire.Fire(fixer)