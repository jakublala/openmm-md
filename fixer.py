from pdbfixer import PDBFixer
from openmm.app import PDBFile
import fire

def main(
        filename=None
):
    if filename is None:
        raise ValueError('Filename is required')
    
    # this script added the C-terminus oxygen atom
    fixer = PDBFixer(filename='S1_Best_A.pdb')
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()

    with open(f'{filename}_fixed.pdb', 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


if __name__ == '__main__':
    fire.Fire(main)