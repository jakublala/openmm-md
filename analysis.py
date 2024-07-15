import mdtraj as md
import fire

def main(
        filename=None
):
    if filename is None:
        raise ValueError('Filename is required')
    

    # Load the PDB file
    traj = md.load_pdb(f'{filename}_traj.pdb')

    # Select non-water atoms
    non_water_ions = traj.topology.select('not water and not (name Na or name Cl)')
    
    # Slice the trajectory to include only non-water atoms
    traj_no_water = traj.atom_slice(non_water_ions)


    # Save the modified PDB file
    traj_no_water.save_pdb(f'{filename}_traj.pdb')


if __name__ == '__main__':
    fire.Fire(main)