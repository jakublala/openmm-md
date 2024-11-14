import os
from src.plumed.cv import get_contact_map
import mdtraj as md
from typing import Literal

import logging
logger = logging.getLogger(__name__)

def assert_config(config):
    if config is None:
        raise ValueError('Config is required')
    if config['type'] is None:
        raise ValueError('Type of PLUMED simulation is required')
    if config['pace'] is None:
        raise ValueError('Pace is required')
    if config['barrier'] is None:
        raise ValueError('Barrier is required')
    if config['temperature'] is None:
        raise ValueError('Temperature is required')

def create_plumed_input(
        filepath, 
        output_dir=None,
        config=None,
        mode: Literal['single-chain', 'two-chain'] = 'single-chain',
        ):
    # only works for two specific CVs!!!
    filename = os.path.basename(filepath).split('.')[0]

    assert_config(config)

    if output_dir is None:
        raise ValueError('Output directory is required')
    
    # get atom ids of atoms in each chain
    topology = md.load(f"{output_dir}/{filename}_fixed.pdb").topology
    
    # TODO: add assertions that the wall is at a lower distance than half of the diagonal

    contact_map = get_contact_map(
        filename=filename, 
        output_dir=output_dir,
        spot1_residues=config['spot1_residues'],
        spot2_residues=config['spot2_residues']
        )
    
    print(contact_map)

    assert 0 == 1


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
        residue_id = topology.atom(atom_id).residue.index + 1
        binder_residues.append(residue_id)
    binder_residues = list(set(binder_residues))
    com_residues_binder = ','.join(f"@CA-A_{i}" for i in binder_residues)

    plumed_content = get_plumed_content(config, output_dir, filename, binder_residues, com_residues_binder, com_residues_binding_site)

    with open(f'{output_dir}/{filename}_plumed.dat', 'w') as f:
        f.write(plumed_content)

    write_pymol_commands(filename, binding_site_residues, contact_residues, output_dir)


def get_plumed_content(
        config, 
        output_dir,
        filename,
        mode: Literal['single-chain', 'two-chain'] = 'single-chain',
        
        ):
    # TODO: rewrite this in a proper PLUMED SCHEMA

    if mode == 'single-chain':
        whole_molecules_content = "WHOLEMOLECULES ENTITY0=@protein"
    elif mode == 'two-chain':
        whole_molecules_content = f"""chain_A: GROUP ATOMS=@protein-A
chain_B: GROUP ATOMS=@protein-B
WHOLEMOLECULES ENTITY0=chain_A ENTITY1=chain_B"""
# this might not work later on when I do two chains

    if mode == 'single-chain':
        com_content = f"""c1: COM ATOMS={com_residues_binder}
c2: COM ATOMS={com_residues_binding_site}
"""
    elif mode == 'two-chain':
        com_content = f"""c1: COM ATOMS=@protein-A
c2: COM ATOMS=@protein-B
"""
        
    

    plumed_content = f"""MOLINFO STRUCTURE={output_dir}/{filename}_fixed.pdb
{whole_molecules_content}
{com_content}
d: DISTANCE ATOMS=c1,c2
cmap: CONTACTMAP ... 
{contact_content}
\tSWITCH={{RATIONAL R_0={config['cutoff']}}}
\tSUM
...
"""
    
    
    if config['type'] == 'opes-explore':
        plumed_content += "opes: OPES_METAD_EXPLORE ...\n"
    elif config['type'] == 'opes':
        plumed_content += "opes: OPES_METAD ...\n"
    else:
        raise ValueError(f"Invalid type: {config['type']}")
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
    return plumed_content


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
