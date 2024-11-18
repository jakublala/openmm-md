from openmm.unit import nanoseconds, picoseconds
def get_checkpoint_interval(timestep):
    return int((1 * nanoseconds) / (timestep * 0.001 * picoseconds)) 


import MDAnalysis as mda
def get_atom_ids_from_chain(chain_id: str, filename: str, output_dir: str):
    return mda.Universe(f'{output_dir}/{filename}_fixed.pdb').select_atoms(f'chainid {chain_id}').ids