import os
from Bio.PDB import PDBParser

import logging
logger = logging.getLogger(__name__)

def create_opes_input(
        filepath, 
        config=None,
        type='opes',
        output_dir=None
        ):
    filename = os.path.basename(filepath).split('.')[0]

    if config is None:
        raise ValueError('Config is required')
    
    if output_dir is None:
        raise ValueError('Output directory is required')

    # get atom ids of atoms in each chain
    import mdtraj as md
    traj = md.load(f"{output_dir}/{filename}_fixed.pdb")
    chain_A_indices = traj.topology.select(f'chainid 0')
    chain_B_indices = traj.topology.select(f'chainid 1')

    # TODO: add assertions that the wall is at a lower distance
    # than half of the diagonal

    # Create atom mapping and get contact indices
    from src.plumed.cv import get_interface_contact_indices
    
    contact_residues = get_interface_contact_indices(filename, cutoff=config['cutoff'], output_dir=output_dir)
    
    # Generate PLUMED input content
    contact_str = '\n'.join(f"\tATOMS{i+1}=@CA-A_{a},@CA-B_{b}" for i, (a, b) in enumerate(contact_residues))

    # assume
    # chain A is the binder
    # chain B is the target

    # assert that chain A is shorter
    assert len(chain_A_indices) < len(chain_B_indices), "Chain A is not the binder"

    binding_site_residues = []
    for _, j in contact_residues:
        # assuming second chain is the target and thus has the binding site
        binding_site_residues.append(j)
    com_residues = [f"@CA-B_{i}" for i in binding_site_residues]
    logger.info(f"The binding site has {len(com_residues)} residues")
    com_residues_binding_site = ','.join(com_residues)


    binder_residues = []
    for i in chain_A_indices:
        atom_id = i
        residue_id = traj.topology.atom(atom_id).residue.index + 1
        binder_residues.append(residue_id)
    binder_residues = list(set(binder_residues))
    com_residues_binder = ','.join(f"@CA-A_{i}" for i in binder_residues)
    
    plumed_content = f"""MOLINFO STRUCTURE={output_dir}/{filename}_fixed.pdb
chain_A: GROUP ATOMS={chain_A_indices[0]+1}-{chain_A_indices[-1]+1}
chain_B: GROUP ATOMS={chain_B_indices[0]+2}-{chain_B_indices[-1]+2}
WHOLEMOLECULES ENTITY0=chain_A ENTITY1=chain_B
c1: COM ATOMS={com_residues_binder}
c2: COM ATOMS={com_residues_binding_site}
d: DISTANCE ATOMS=c1,c2
cmap: CONTACTMAP ... 
{contact_str}
\tSWITCH={{RATIONAL R_0={config['cutoff']}}}
\tSUM
...
"""
    if type == 'opes_explore':
        plumed_content += "opes: OPES_METAD_EXPLORE ...\n"
    elif type == 'opes':
        plumed_content += "opes: OPES_METAD ...\n"
    else:
        raise ValueError(f"Invalid type: {type}")
    plumed_content += f"""\tARG=cmap,d PACE={config['pace']} BARRIER={config['barrier']}
\tTEMP={config['temperature']}
\tFILE={output_dir}/{filename}.kernels
"""
    if config['restart_rfile'] is not None:
        plumed_content += f"\tSTATE_RFILE={config['restart_rfile']}\n\tRESTART=YES\n"
    else:
        pass

    plumed_content += f"""\tSTATE_WFILE={output_dir}/{filename}.state\n\tSTATE_WSTRIDE={config['state_wstride']}
...
uwall: UPPER_WALLS ARG=d AT={config['upper_wall.at']} KAPPA=150.0 EXP={config['upper_wall.exp']} EPS=1 OFFSET=0
PRINT ARG=cmap,d,opes.*,uwall.bias STRIDE={config['stride']} FILE={output_dir}/{filename}.colvar
"""
    with open(f'{output_dir}/{filename}_plumed.dat', 'w') as f:
        f.write(plumed_content)

    write_pymol_commands(filename, binding_site_residues, contact_residues, output_dir)




def write_pymol_commands(filename, binding_site_residues, contact_residues, output_dir):
    # Write PyMOL visualization commands
    pymol_commands = f"""load {output_dir}/{filename}_fixed.pdb

# Color chain B (target) in gray first
select target, chain B
color gray50, target

# Then color binding site residues yellow
select binding_site, resi {'+'.join(str(i) for i in binding_site_residues)} and chain B
color green, binding_site

# Color chain A (binder) in blue
select binder, chain A
color marine, binder

# Color contact residues in chain A in red
select binder_contacts, resi {'+'.join(str(i) for i, _ in contact_residues)} and chain A
color red, binder_contacts

# Set nice visualization
hide everything
show cartoon
"""

    with open(f'{output_dir}/{filename}_pymol_commands.txt', 'w') as f:
        f.write(pymol_commands)
