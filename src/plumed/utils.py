import pandas as pd
from openmm.unit import nanoseconds, picoseconds
def get_checkpoint_interval(timestep):
    return int((1 * nanoseconds) / (timestep * 0.001 * picoseconds)) 


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

def get_last_checkpoint_timestep(out_file, checkpoint_interval):
    df_out = read_openmm_output(out_file)
    last_timestep = df_out['Step'].iloc[-1]
    num_checkpoints = last_timestep // checkpoint_interval
    return num_checkpoints * checkpoint_interval

def process_hills_for_restart(hills_file, time_of_last_hill):
    filtered_lines = []
    with open(hills_file, 'r') as f:
        # Keep header lines (starting with #) and hills up to time_of_last_hill
        for line in f:
            if line.startswith('#'):
                filtered_lines.append(line)
                continue
            
            # Hills file columns are space-separated, time is in the first column
            try:
                time = float(line.split()[0])
                if time <= time_of_last_hill:
                    filtered_lines.append(line)
            except (ValueError, IndexError):
                continue
    
    return filtered_lines



def assert_correct_metad_grid(config):
    grid_min = [int(x) for x in config['metad.grid_min'].split(',')]
    grid_max = [int(x) for x in config['metad.grid_max'].split(',')]
    grid_bin = [int(x) for x in config['metad.grid_bin'].split(',')]
    sigmas = [float(x) for x in config['metad.sigma'].split(',')]

    for i, _ in enumerate(config['cvs']):
        grid_spacing = (grid_max[i] - grid_min[i]) / grid_bin[i]
        assert grid_spacing / 2 < sigmas[i], f'Grid spacing for {config['cvs'][i]} is too large, increase SIGMA'

        