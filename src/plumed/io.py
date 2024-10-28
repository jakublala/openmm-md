import os
from Bio.PDB import PDBParser

def create_opes_input(
        filepath, 
        cv_string, 
        config=None,
        type='opes'
        ):
    filename = os.path.basename(filepath).split('.')[0]

    parser = PDBParser()
    structure = parser.get_structure('protein', f'tmp/{filename}/{filename}_fixed.pdb')

    # get the initial atom ID for each chain
    # get the IDs for chain A and B
    atom_ids_A = [atom.get_serial_number() for atom in structure[0]['A'].get_atoms()]
    atom_ids_B = [atom.get_serial_number() for atom in structure[0]['B'].get_atoms()]

    if config is None:
        raise ValueError('Config is required')

    # get atoms that are part of the binding site on the target protein
    # and have been included as the part of the contact map
    from src.plumed.cv import get_interface_contact_indices
    contact_indices = get_interface_contact_indices(filename, cutoff=config['cutoff'])

    print(contact_indices)
    assert 0 == 1


    with open(f'tmp/{filename}/{filename}_plumed.dat', 'w') as f:
        content = f"""MOLINFO STRUCTURE=tmp/{filename}/{filename}_fixed.pdb
chain_A: GROUP ATOMS={atom_ids_A[0]}-{atom_ids_A[-1]}
chain_B: GROUP ATOMS={atom_ids_B[0]}-{atom_ids_B[-1]}
WHOLEMOLECULES ENTITY0=chain_A ENTITY1=chain_B
c1: COM ATOMS=chain_A
c2: COM ATOMS=chain_B
d: DISTANCE ATOMS=c1,c2
cmap: CONTACTMAP ... 
{cv_string}
\tSWITCH={{RATIONAL R_0={config['cutoff']}}}
\tSUM
...
"""
        if type == 'opes_explore':
            content += "opes: OPES_METAD_EXPLORE ...\n"
        elif type == 'opes':
            content += "opes: OPES_METAD ...\n"
        else:
            raise ValueError(f"Invalid type: {type}")
        content += f"""\tARG=cmap,d PACE={config['pace']} BARRIER={config['barrier']}
\tTEMP={config['temperature']}
\tFILE=tmp/{filename}/{filename}.kernels
"""
        if config['restart_rfile'] is not None:
            content += f"\tSTATE_RFILE={config['restart_rfile']}\n\tRESTART=YES\n"
        else:
            pass

        content += f"""\tSTATE_WFILE=tmp/{filename}/{filename}.state\n\tSTATE_WSTRIDE={config['chk_stride']}
...
PRINT ARG=cmap,d,opes.* STRIDE={config['stride']} FILE=tmp/{filename}/{filename}.colvar
"""

        f.write(content)

