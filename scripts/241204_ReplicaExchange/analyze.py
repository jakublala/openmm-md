from openmmtools.multistate import ParallelTemperingAnalyzer
import pathlib

def main():
    storage_path = pathlib.Path('results/replica_exchange.nc')
    import netCDF4
    nc = netCDF4.Dataset(storage_path, 'r')
    print(nc.variables.keys())
    n_replicas = nc['states'].shape[0]
    n_iterations = nc['states'].shape[1]

    from openmm.app import PDBFile
    from openmm.unit import nanometer, angstrom
    pdb_path = '../../data/241010_FoldingUponBinding/output/SUMO-1C/241128-MetaD/sumo1c_fixed.pdb'
    # copy pdb file to here directory
    import shutil
    shutil.copy(pdb_path, 'results/replica_0.pdb')
    
    # loading traj and write to dcd format
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader
    positions = nc['positions'][:].swapaxes(0, 1)
    positions = positions * 10  # Convert nm to Angstroms for MDAnalysis

    # Extract and convert box vectors (assuming they're in nm)
    box_vectors = nc['box_vectors'][:] * 10  # Convert nm to Angstroms
    
    # Create trajectory coordinates in the format MDAnalysis expects
    coordinates = positions[0]  # Taking first replica
    trajectory = coordinates.reshape(len(coordinates), -1, 3)
    
    # Create universe with the correct trajectory format
    u = mda.Universe('results/replica_0.pdb', trajectory, format=MemoryReader)

    # Set dimensions for each frame
    for ts, box in zip(u.trajectory, box_vectors[0]):
        # Convert box vectors to lengths and angles
        a = np.linalg.norm(box[0])
        b = np.linalg.norm(box[1])
        c = np.linalg.norm(box[2])
        
        # Set box dimensions
        alpha = np.degrees(np.arccos(np.dot(box[1], box[2]) / (b * c)))
        beta = np.degrees(np.arccos(np.dot(box[0], box[2]) / (a * c)))
        gamma = np.degrees(np.arccos(np.dot(box[0], box[1]) / (a * b)))
        ts.dimensions = [a, b, c, alpha, beta, gamma]

    # Write the centered trajectory
    with mda.Writer("results/replica_0.dcd", u.select_atoms("protein").n_atoms) as W:
        for ts in u.trajectory:
            W.write(u.select_atoms("protein"))


    # plot first atom from positions with matplitlbi
    import matplotlib.pyplot as plt
    plt.plot(positions[:, 0, 0])
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.savefig('results/first_atom_positions.png')

    # plot velocities
    velocities = nc['velocities'][:].swapaxes(0, 1)
    plt.plot(velocities[:, 0, 0])
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.savefig('results/first_atom_velocities.png')


    print(f"{nc['states'][:]=}")
    print(f"{nc['accepted'][:]=}")
    print(f"{nc['proposed'][:]=}")
    
    # plot the states
    plt.plot(nc['states'][:])
    plt.xlabel('Time')
    plt.ylabel('State')
    plt.savefig('results/states.png')

    from openmmtools.multistate import MultiStateReporter
    reporter = MultiStateReporter('results/replica_exchange.nc', open_mode='r')
    print(f"{reporter.read_sampler_states(0)[0].positions=}")

    assert 0 == 1
    # plot the energies
    print(f"{nc['energies'][:]}")
    plt.plot(nc['energies'][:])
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.savefig('results/energies.png')


if __name__ == "__main__":
    main()