import pandas as pd
from openmm.unit import nanoseconds, picoseconds

import MDAnalysis as mda
def get_atom_ids_from_chain(chain_id: str, filename: str, output_dir: str):
    return mda.Universe(f'{output_dir}/{filename}_fixed.pdb').select_atoms(f'chainid {chain_id}').ids


def get_pace_from_metad(plumed_file):
    with open(plumed_file, 'r') as f:
        for line in f:
            if 'PACE=' in line:
                # Find PACE value using string split
                pace_str = line.split('PACE=')[1].split()[0].strip(',')
                return int(pace_str)
    return None  # Return None if PACE not found

def read_openmm_output(out_file: str) -> pd.DataFrame:
    # Read the file, skip the first line (header starts with #)
    df = pd.read_csv(out_file, 
                     sep='\t',
                     comment='#',
                     names=['Progress', 'Step', 'Time', 'Potential_Energy', 
                           'Kinetic_Energy', 'Total_Energy', 'Temperature',
                           'Box_Volume', 'Density', 'Speed', 'Time_Remaining'])
    
    # Clean up the columns
    df['Progress'] = df['Progress'].str.rstrip('%').astype(float) / 100
    df['Time'] = df['Time'].astype(float)  # Time in ps
    
    return df

def assert_correct_metad_grid(config):
    # CV1
    assert config['cv1.type'] is not None, "Expected CV1 type"
    assert config['cv1.grid_min'] is not None, "Expected CV1 grid min"
    assert config['cv1.grid_max'] is not None, "Expected CV1 grid max"
    assert config['cv1.grid_bin'] is not None, "Expected CV1 grid bin"
    assert config['cv1.sigma'] is not None, "Expected CV1 sigma"

    grid_spacing = (config['cv1.grid_max'] - config['cv1.grid_min']) / config['cv1.grid_bin']
    assert grid_spacing / 2 < config['cv1.sigma'], f'Grid spacing for {config['cv1.type']} is too large, increase SIGMA'
    
    # CV2
    assert config['cv2.type'] is not None, "Expected CV2 type"
    assert config['cv2.grid_min'] is not None, "Expected CV2 grid min"
    assert config['cv2.grid_max'] is not None, "Expected CV2 grid max"
    assert config['cv2.grid_bin'] is not None, "Expected CV2 grid bin"
    assert config['cv2.sigma'] is not None, "Expected CV2 sigma"

    grid_spacing = (config['cv2.grid_max'] - config['cv2.grid_min']) / config['cv2.grid_bin']
    assert grid_spacing / 2 < config['cv2.sigma'], f'Grid spacing for {config['cv2.type']} is too large, increase SIGMA'



def prepare_plumed_file_for_restart(plumed_file, output_dir, system):
    with open(plumed_file, 'r') as f:
        lines = f.readlines()
    
    # change the RESTART=NO to RESTART=YES
    for i, line in enumerate(lines):
        if 'RESTART=NO' in line:
            lines[i] = line.replace('RESTART=NO', 'RESTART=YES')
            break

    # find FILE=<anything>.hills and FILE=<anything>.colvar and replace them with FILE=<output_dir>/restart-<index>/<anything>.hills and FILE=<output_dir>/restart-<index>/<anything>.colvar
    # Use regex to find and replace FILE= patterns for both .hills and .colvar files
    import re
    for i, line in enumerate(lines):
        if 'FILE=' in line:
            # Replace any FILE=<path>/<name>.hills or FILE=<path>/<name>.colvar
            line = re.sub(r'FILE=.*?\.hills', f'FILE={output_dir}/{system}.hills', line)
            line = re.sub(r'FILE=.*?\.colvar', f'FILE={output_dir}/{system}.colvar', line)
            lines[i] = line

    with open(f'{output_dir}/{system}_plumed.dat', 'w') as f:
        f.writelines(lines)

    return f'{output_dir}/{system}_plumed.dat'