from openmmtools.multistate import ParallelTemperingAnalyzer
import pathlib
import shutil
import numpy as np
from openmm.app import PDBFile
from openmm.unit import nanometer, angstrom
from MDAnalysis.coordinates.memory import MemoryReader
import MDAnalysis as mda
from tqdm import tqdm

def process_trajectory(nc, pdb_path, N_REPLICAS):
    for i in tqdm(range(N_REPLICAS), desc="Producing DCD trajectories for replicas"):
        dcd_filepath = f'results/replica_{i}.dcd'
        
        shutil.copy(pdb_path, f'results/replica_{i}.pdb')
        positions = nc['positions'][:].swapaxes(0, 1)
        positions = positions * 10  # Convert nm to Angstroms for MDAnalysis

        # Extract and convert box vectors (assuming they're in nm)
        box_vectors = nc['box_vectors'][:] * 10  # Convert nm to Angstroms
        
        # Create trajectory coordinates in the format MDAnalysis expects
        coordinates = positions[0]  # Taking first replica
        trajectory = coordinates.reshape(len(coordinates), -1, 3)
        
        # Create universe with the correct trajectory format
        u = mda.Universe('results/replica_0.pdb', trajectory, format=MemoryReader, guess_bonds=True)

        # Calculate box dimensions once since they're constant
        box = box_vectors[0, 0]  # Take first frame of first replica
        a = np.linalg.norm(box[0])
        b = np.linalg.norm(box[1])
        c = np.linalg.norm(box[2])
        alpha = np.degrees(np.arccos(np.dot(box[1], box[2]) / (b * c)))
        beta = np.degrees(np.arccos(np.dot(box[0], box[2]) / (a * c)))
        gamma = np.degrees(np.arccos(np.dot(box[0], box[1]) / (a * b)))
        dimensions = [a, b, c, alpha, beta, gamma]
            
        # Set dimensions for each timestep
        for ts in u.trajectory:
            ts.dimensions = dimensions
        
        protein = u.select_atoms('protein')
        
        # Use MDAnalysis transformations in the correct order
        from MDAnalysis.transformations import wrap, unwrap, center_in_box
        
        workflow = [
            unwrap(protein), 
            center_in_box(protein, center='geometry'), 
            wrap(protein, compound='segments'),  
            center_in_box(protein, center='mass'),
        ]
        
        u.trajectory.add_transformations(*workflow)

        # Write the centered trajectory
        with mda.Writer(dcd_filepath, u.select_atoms("protein").n_atoms) as W:
            for ts in u.trajectory:
                W.write(u.select_atoms("protein"))
        

def print_acceptance_statistics(nc, n_replicas):
    possible_swaps = {
        (i, i+1): 0 for i in range(n_replicas - 1)
    }

    num_attempts = np.sum(nc['proposed'][1:], axis=0)
    num_accepted = np.sum(nc['accepted'][1:], axis=0)
    acceptance_rate = num_accepted / num_attempts
    for i, j in possible_swaps.keys():
        print(f"Swap {i} <-> {j} Acceptance Rate: {100 * acceptance_rate[i, j]:2f}%")
        
import netCDF4
def main():
    storage_path = pathlib.Path('results/replica_exchange.nc')
    
    nc = netCDF4.Dataset(storage_path, 'r')
    print(nc.variables.keys())
    N_REPLICAS = nc['states'].shape[1]
    n_iterations = nc['states'].shape[1]

    pdb_path = '../../data/241010_FoldingUponBinding/output/SUMO-1C/241128-MetaD/sumo1c_fixed.pdb'
    # process_trajectory(nc, pdb_path, N_REPLICAS)
    # compute the statistics of accepted swaps

    print_acceptance_statistics(nc, N_REPLICAS)

    # # plot first atom from positions with matplitlbi
    # import matplotlib.pyplot as plt
    # plt.plot(positions[:, 0, 0])
    # plt.xlabel('Time')
    # plt.ylabel('Position')
    # plt.savefig('results/first_atom_positions.png')

    # # plot velocities
    # velocities = nc['velocities'][:].swapaxes(0, 1)
    # plt.plot(velocities[:, 0, 0])
    # plt.xlabel('Time')
    # plt.ylabel('Velocity')
    # plt.savefig('results/first_atom_velocities.png')


    # print(f"{nc['states'][:]=}")
    # print(f"{nc['accepted'][:]=}")
    # print(f"{nc['proposed'][:]=}")
    
    # # plot the states
    # plt.plot(nc['states'][:])
    # plt.xlabel('Time')
    # plt.ylabel('State')
    # plt.savefig('results/states.png')

    # # from openmmtools.multistate import MultiStateReporter
    # # reporter = MultiStateReporter('results/replica_exchange.nc', open_mode='r')
    # # print(f"{reporter.read_sampler_states(0)[0].positions=}")

    
    # # plot the energies
    # print(f"{nc['energies'][:]}")
    # plt.plot(nc['energies'][:])
    # plt.xlabel('Time')
    # plt.ylabel('Energy')
    # plt.savefig('results/energies.png')


if __name__ == "__main__":
    main()