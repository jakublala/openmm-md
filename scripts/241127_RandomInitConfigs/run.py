import re
import fire
import os
import shutil
import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations import unwrap, center_in_box

from src.analysis.utils import get_file_by_extension
from src.analysis.io import read_colvar_file

TIMESTEP = 0.002 # in ps
STRIDE = 500
LOGGING_FREQUENCY = 100 # of STRIDE

def main(system: str, date: str, project: str, num_configs: int = 5):
    plumed_file = get_file_by_extension(f"../../data/{project}/output/{system}/{date}", '.dat')
    dcd_file = get_file_by_extension(f"../../data/{project}/output/{system}/{date}", '.dcd')
    pdb_file = get_file_by_extension(f"../../data/{project}/output/{system}/{date}", 'fixed.pdb')
    colvar_file = get_file_by_extension(f"../../data/{project}/output/{system}/{date}", '.colvar')
    colvar_df = read_colvar_file(colvar_file)
    # find min and max value of 'd' colvar
    min_d = colvar_df['d'].min()
    max_d = colvar_df['d'].max()
    
    d_values = np.linspace(min_d, max_d, num_configs)
    
    universe = mda.Universe(pdb_file, dcd_file)
    # Add both unwrapping and centering transformations
    universe.trajectory.add_transformations(
        unwrap(universe.atoms),
        center_in_box(universe.select_atoms('chainid B'), center='geometry', wrap=False)
    )

    # get unit cell dimension
    unit_cell_dimensions = universe.atoms.dimensions
    assert unit_cell_dimensions[0] == unit_cell_dimensions[1] == unit_cell_dimensions[2], "Box is not cubic"
    d_box = unit_cell_dimensions[0]


    # make a directory for the new configs
    directory = f"../../data/{project}/output/{system}/{date}/extracted_frames_from_traj"
    os.makedirs(directory, exist_ok=True)

    colvar_traj_values = []
    # find d values from colvar_df that match the trajectory frames
    for i in range(universe.trajectory.n_frames):
        colvar_traj_values.append(colvar_df.iloc[i * LOGGING_FREQUENCY].d)
    
    closest_frame_ids = []
    for d in d_values:
        # stride was 500 in the original simulation
        # timestep was 2 fs in the original simulation
        # we got a colvar frame, every 1 ps
        # we got a frame from the trajectory every 100 * 500 * timesteps
        # traj_dt = 100 * 500 * TIMESTEP # for this case, 100 ps

        # find the index of the closest d value to the given d
        closest_frame_ids.append(np.abs(np.array(colvar_traj_values) - d).argmin())

    

    for i, frame_id in enumerate(closest_frame_ids):
        output_filename = f"{directory}/CD28-G-{round(d_values[i], 1)}.pdb"
        # replace the closest frame with the original fixed.pdb file
        if i == 0:
            # just copy the fixed.pdb file over
            shutil.copy(pdb_file, output_filename)
            continue


        # Set the current frame in the universe
        universe.trajectory[frame_id]
        
        chain_A_com = universe.select_atoms('chainid A').center_of_mass()
        chain_B_com = universe.select_atoms('chainid B').center_of_mass()
        v_com = chain_A_com - chain_B_com
        translate_v = np.round(v_com / d_box, 0) * d_box # for chain B

        # translate all chain B atoms in this timestep
        universe.select_atoms('chainid B').translate(translate_v)
        
        # Write the PDB for the current frame
        with mda.Writer(output_filename) as pdb:
            pdb.write(universe.atoms)


        assert_consistent_colvar_d(universe, d_values[i], plumed_file)

from collections import OrderedDict
def get_coms_from_plumed_file(plumed_file):
    # TODO: make it work with single-chain systems
    com_groups = OrderedDict({'A': [], 'B': []})
    
    with open(plumed_file, 'r') as file:
        for line in file:
            if not line.startswith(('c1:', 'c2:')):
                continue
                
            residues = re.findall(r'@CA-([AB])_(\d+)', line)
            for chain, num in residues:
                com_groups[chain].append(int(num))
    
    return {chain: sorted(nums) for chain, nums in com_groups.items()}

def assert_consistent_colvar_d(universe, d, plumed_file):
    com_groups = get_coms_from_plumed_file(plumed_file)
    print(com_groups)
    # Use spaces instead of commas to separate residue IDs
    chain_A_com = universe.select_atoms(f'name CA and chainid A and resid {" ".join(map(str, com_groups["A"]))}').center_of_mass()
    chain_B_com = universe.select_atoms(f'name CA and chainid B and resid {" ".join(map(str, com_groups["B"]))}').center_of_mass()

    print(universe.select_atoms(f'name CA and chainid A and resid {" ".join(map(str, com_groups["A"]))}').positions)

    print(chain_A_com, chain_B_com)

    d_COM = np.linalg.norm(chain_A_com - chain_B_com)
    assert np.isclose(d_COM / 10, d), f"d_COM is not consistent with d: {d_COM / 10} != {d}"

if __name__ == "__main__":
    fire.Fire(main)