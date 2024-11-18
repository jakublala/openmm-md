import os
from src.plumed.cv import get_contact_map
from typing import Literal

from src.models import ContactMap, Segment
from typing import Optional

from src.plumed.utils import get_atom_ids_from_chain

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
    
    # TODO: add assertions that the wall is at a lower distance than half of the diagonal

    contact_map = get_contact_map(
        filename=filename, 
        output_dir=output_dir,
        spot1_residues=config['spot1_residues'],
        spot2_residues=config['spot2_residues'],
        mode=mode
        )
    
    plumed_content = get_plumed_content(
        config=config, 
        output_dir=output_dir, 
        filename=filename, 
        contact_map=contact_map,
        mode=mode,
        spot1_segment=config['spot1_residues'],
        spot2_segment=config['spot2_residues']
        )

    with open(f'{output_dir}/{filename}_plumed.dat', 'w') as f:
        f.write(plumed_content)

    write_pymol_commands(filename, output_dir, contact_map)


def get_plumed_content(
        config: dict, 
        output_dir: str,
        filename: str,
        contact_map: ContactMap,
        mode: Literal['single-chain', 'two-chain'] = 'single-chain',
        spot1_segment: Optional[Segment] = None,
        spot2_segment: Optional[Segment] = None,
        ) -> str:
    # TODO: rewrite this in a proper PLUMED SCHEMA

    if spot1_segment is not None and spot2_segment is not None:
        assert spot1_segment.indexing == spot2_segment.indexing == 1, (
            f"Expected indexing 1 for both spot1 and spot2, got {spot1_segment.indexing} and {spot2_segment.indexing}"
        )
    spot1_com_CAs = sorted(set([str(c.residue1) for c in contact_map.contacts]))
    spot2_com_CAs = sorted(set([str(c.residue2) for c in contact_map.contacts]))
    spot1_com_CAs = ",".join(spot1_com_CAs)
    spot2_com_CAs = ",".join(spot2_com_CAs)

    if mode == 'single-chain':
        whole_molecules_content = "WHOLEMOLECULES ENTITY0=@protein"
    elif mode == 'two-chain':
        atom_ids_chain_A = get_atom_ids_from_chain('A', filename, output_dir)
        atom_ids_chain_B = get_atom_ids_from_chain('B', filename, output_dir)
        group1_content = f"chain_A: GROUP ATOMS={",".join([f"{i}" for i in atom_ids_chain_A])}"
        group2_content = f"chain_B: GROUP ATOMS={",".join([f"{i}" for i in atom_ids_chain_B])}"

        whole_molecules_content = f"""{group1_content}
{group2_content}
WHOLEMOLECULES ENTITY0=chain_A ENTITY1=chain_B"""
    com_content = f"""c1: COM ATOMS={spot1_com_CAs}
c2: COM ATOMS={spot2_com_CAs}"""
    
    plumed_content = f"""MOLINFO STRUCTURE={output_dir}/{filename}_fixed.pdb
{whole_molecules_content}
{com_content}
d: DISTANCE ATOMS=c1,c2
cmap: CONTACTMAP ... 
{contact_map}
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


def write_pymol_commands(filename, output_dir, contact_map: ContactMap):
    # assuming indexing of 1
    # TODO: do proper assertion
    residues_site_1 = set([i.residue1.index for i in contact_map.contacts])
    residues_site_2 = set([i.residue2.index for i in contact_map.contacts])
    binding_site_1_str = '+'.join(str(i) for i in residues_site_1)
    binding_site_2_str = '+'.join(str(i) for i in residues_site_2)
    site_1_chain_id = contact_map.contacts[0].residue1.chain_id
    site_2_chain_id = contact_map.contacts[0].residue2.chain_id
    # Write PyMOL visualization commands
    pymol_commands = f"""load {output_dir}/{filename}_fixed.pdb

# First color everything blue
color blue, all

# Then color binding site residues yellow
select binding_site, resi {binding_site_1_str} and chain {site_1_chain_id}
color green, binding_site

# Color contact residues in chain A in red
select binder_contacts, resi {binding_site_2_str} and chain {site_2_chain_id}
color red, binder_contacts

# Set nice visualization
hide everything
show cartoon
"""

    with open(f'{output_dir}/{filename}_pymol_commands.txt', 'w') as f:
        f.write(pymol_commands)
