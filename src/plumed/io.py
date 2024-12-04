import os
from src.plumed.cv import get_contact_map
from typing import Literal
import numpy as np

from src.models import ContactMap, Segment
from typing import Optional

from src.analysis.utils import get_file_by_extension
from src.plumed.utils import get_atom_ids_from_chain

import logging
logger = logging.getLogger(__name__)

def assert_config(config):
    if config is None:
        raise ValueError('Config is required')
    if config['type'] is None:
        raise ValueError('Type of PLUMED simulation is required')
    if config['temperature'] is None:
        raise ValueError('Temperature is required')

def create_plumed_input(
        filename, 
        output_dir=None,
        config=None,
        mode: Literal['single-chain', 'two-chain'] = 'single-chain',
        ):
    # only works for two specific CVs!!!
    
    assert_config(config)

    get_file_by_extension(output_dir, 'equilibrated.pdb', assert_exists=True)
    
    if output_dir is None:
        raise ValueError('Output directory is required')
    
    # TODO: add assertions that the wall is at a lower distance than half of the diagonal

    contact_map = get_contact_map(
        filename=filename, 
        cutoff=config['cutoff'],
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
        range_atom_ids_chain_A = [atom_ids_chain_A[0], atom_ids_chain_A[-1]]
        range_atom_ids_chain_B = [atom_ids_chain_B[0], atom_ids_chain_B[-1]]
        assert np.all(np.arange(range_atom_ids_chain_A[0], range_atom_ids_chain_A[1] + 1) == atom_ids_chain_A), "Atom ids in chain A are not continuous"
        assert np.all(np.arange(range_atom_ids_chain_B[0], range_atom_ids_chain_B[1] + 1) == atom_ids_chain_B), "Atom ids in chain B are not continuous"
        group1_content = f"chain_A: GROUP ATOMS={range_atom_ids_chain_A[0]}-{range_atom_ids_chain_A[1]}"
        group2_content = f"chain_B: GROUP ATOMS={range_atom_ids_chain_B[0]}-{range_atom_ids_chain_B[1]}"

    assert len(config['cvs']) == 2, "Expected two CVs, got {len(config['cvs'])}. Don't support any more or less."
    
    def resolve_cv_content(cv: str, config: dict):
        # TODO: make this into proper dataclasses that are resolved with __str__
        if cv == 'd':
            return "d: DISTANCE ATOMS=c1,c2\n"
        elif cv == 'cmap':
            return f"""cmap: CONTACTMAP ... 
{contact_map}
\tSWITCH={{RATIONAL R_0={config['cutoff']}}}
\tSUM
...
"""
        elif cv == 'sasa':
            sasa_algo = config['sasa.algo'].upper()
            assert sasa_algo in ['HASEL', 'LCPO'], f"Invalid SASA algorithm: {sasa_algo}. Expected 'HASEL' or 'LCPO'."
            assert config['sasa.spot_id'] in [1, 2], f"Invalid SASA spot id: {config['sasa.spot_id']}. Expected 1 or 2."
            if config['sasa.spot_id'] == 1:
                residues = contact_map.all_residues1
            else:
                residues = contact_map.all_residues2
            atom_ids = sorted(set([atom_id for res in residues for atom_id in res.atom_indices]))
            atom_ids_str = ",".join(str(atom_id) for atom_id in atom_ids)
            logger.info(f"SASA in total made up from {len(atom_ids)} atoms")
            return f"sasa: SASA_{sasa_algo} TYPE=TOTAL ATOMS={atom_ids_str} NL_STRIDE={config['sasa.stride']}"
        else:
            raise NotImplementedError(f"Unsupported CV: {cv}")

    cv1_content = resolve_cv_content(config['cvs'][0], config)
    cv2_content = resolve_cv_content(config['cvs'][1], config)
    
    # TODO: I don't think this is correct for a single-chain
    whole_molecules_content = f"""{group1_content}
{group2_content}
WHOLEMOLECULES ENTITY0=chain_A ENTITY1=chain_B"""
    com_content = f"""c1: COM ATOMS={spot1_com_CAs}
c2: COM ATOMS={spot2_com_CAs}"""

    plumed_content = f"""MOLINFO STRUCTURE={output_dir}/{filename}_fixed.pdb
{whole_molecules_content}
{com_content}
{cv1_content}
{cv2_content}
"""
    
    cv_arg_content = ','.join(config['cvs'])
    
    
    if config['type'] == 'opes-explore':
        plumed_content += f"""opes: OPES_METAD_EXPLORE ...
\tARG={cv_arg_content} PACE={config['opes.pace']} BARRIER={config['opes.barrier']}
\tTEMP={config['temperature']}
\tFILE={output_dir}/{filename}.kernels
\tSTATE_WFILE={output_dir}/{filename}.state
\tSTATE_WSTRIDE={config['state_wstride']}
...
"""
    elif config['type'] == 'opes':
        plumed_content += f"""opes: OPES_METAD ...
\tARG={cv_arg_content} PACE={config['opes.pace']} BARRIER={config['opes.barrier']}
\tTEMP={config['temperature']}
\tFILE={output_dir}/{filename}.kernels
\tSTATE_WFILE={output_dir}/{filename}.state
\tSTATE_WSTRIDE={config['state_wstride']}
\tSIGMA_MIN=0.001,0.001
...
"""
    elif config['type'] == 'metad':
        # assert that if grid, there's no artifacts
        assert_correct_metad_grid(config)

        restart_str = 'YES' if config['restart'] else 'NO'
        plumed_content += f"""metad: METAD ...
\tARG={cv_arg_content} PACE={config['metad.pace']} SIGMA={config['metad.sigma']} HEIGHT={config['metad.height']}
\tGRID_MIN={config['metad.grid_min']} GRID_MAX={config['metad.grid_max']} GRID_BIN={config['metad.grid_bin']}
\tTEMP={config['temperature']} BIASFACTOR={config['metad.biasfactor']}
\tFILE={output_dir}/{filename}.hills
\tRESTART={restart_str}
...
"""
    else:
        raise ValueError(f"Invalid type: {config['type']}")


    if config['restart_rfile'] is not None:
        plumed_content += f"\tSTATE_RFILE={config['restart_rfile']}\n\tRESTART=YES\n"
    else:
        pass

    if 'spring' in config and config['spring']:
        plumed_content += f"restraint: RESTRAINT ARG=d AT=0 KAPPA={config['spring.k']}\n"

    type_content = 'opes' if 'opes' in config['type'] else 'metad'
    print_arg = f"{cv_arg_content},{type_content}.*,uwall.bias"
    if 'spring' in config and config['spring']:
        print_arg += ",restraint.bias"

    plumed_content += f"""uwall: UPPER_WALLS ARG=d AT={config['upper_wall.at']} KAPPA={config['upper_wall.kappa']} EXP={config['upper_wall.exp']} EPS=1 OFFSET=0
PRINT ARG={print_arg} STRIDE={config['stride']} FILE={output_dir}/{filename}.colvar
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
