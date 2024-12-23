from src.analysis.metad import plot_hills
from src.analysis.colvar import plot_colvar_trajectories

import matplotlib.pyplot as plt

from src.analysis.utils import get_file_by_extension
from src.analysis.io import read_hills_file, read_colvar_file
from src.analysis.colvar import resolve_cv_label

import numpy as np
import netCDF4

TIMESTEP = 0.002 * 10**-3  # in ns
OPES_PACE = 500
TEMP = 300
STRIDE = 500


def plot_replica_analysis(directory, system, cvs, n_replicas=4):
    # Create figure with constrained layout for better spacing
    fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    
    # Add title with padding
    date = directory.split('/')[-1]
    fig.suptitle(f"Replica Analysis for {system}/{date}", fontsize=16, y=1.02)

    # First plot: plot distribution of potential energies
    # load the replica_exchange.nc file with 
    replica_exchange_file = get_file_by_extension(directory, 'replica_exchange.nc')
    nc = netCDF4.Dataset(replica_exchange_file, 'r')

    # print(nc.variables.keys())
    # print(nc.variables['neighborhoods'][:])
    # print(nc.variables['metadata'])
    # print(nc.variables['options'])
    # print(nc.variables['analysis_particle_indices'])



    # float64 energies(iteration, replica, state)
    # TODO: I am not sure what is the difference betweeen the replica and the state dimension
    # HACK: for now, let's just assume that columns are constant-ish
    import pdb; pdb.set_trace()
    potential_energies = np.array(nc.variables['energies'])[:, 0, :].squeeze()
    for i in range(n_replicas):
        axs[0].hist(potential_energies[:, i], bins=100, label=f"Replica {i+1}")
    axs[0].set_xlabel("Potential Energy [reduced units]")
    axs[0].set_ylabel("Frequency")
    axs[0].set_title("Distribution of Potential Energies")
    axs[0].legend()

    # plot the acceptance probability of each type of swap
    possible_swaps = {
            (i, i+1): 0 for i in range(n_replicas - 1)
        }
    
    num_attempts = np.cumsum(nc['proposed'][1:], axis=0)
    num_accepted = np.cumsum(nc['accepted'][1:], axis=0)

    acceptance_rate = num_accepted / num_attempts
    for i, j in possible_swaps.keys():
        axs[1].plot(acceptance_rate[:, i, j], label=f"Swap {i} <-> {j}")
    axs[1].set_xlabel("Iteration")
    axs[1].set_ylabel("Acceptance Rate")
    axs[1].set_title("Acceptance Rate of Each Swap")
    axs[1].legend()

    # plot the states index across the temperatures
    T_MIN = 300
    T_MAX = 600
    TEMPERATURES = [T_MIN + (T_MAX - T_MIN) * (np.exp(float(i) / float(n_replicas-1)) - 1.0) / (np.e - 1.0) for i in range(n_replicas)]
    states = np.array(nc.variables['states'])[:]
    for i in range(states.shape[0]):
        for j in range(states.shape[1]):
            states[i, j] = TEMPERATURES[states[i, j]]
    for i in range(n_replicas):
        axs[2].plot(states[:, i], label=f"Replica {i+1}")
    axs[2].set_xlabel("Iteration")
    axs[2].set_ylabel("Temperature [K]")
    axs[2].set_title("States Index Across Temperatures")
    axs[2].legend()


    # Make all axes square and equal
    for ax in axs.flat:
        ax.set_aspect('auto')
        ax.grid(True, linestyle='--', alpha=0.3)

    # Save with high quality settings
    plt.savefig(
        f"{directory}/{system}_replica_analysis.png",
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    plt.close()    






def plot_replica_trajectory_analysis(directory, system, cvs, n_replicas=4):
    """Create a 3x4 subplot figure showing hills and colvar trajectories for replica analysis.
    
    Args:
        directory (str): Path to the directory containing simulation data
        system (str): Name of the system
        cvs (list): List of collective variables names
    """
    # Create figure with constrained layout for better spacing
    fig, axs = plt.subplots(5, n_replicas, figsize=(16, 20), constrained_layout=True)
    
    # Add title with padding
    date = directory.split('/')[-1]
    fig.suptitle(f"Replica Analysis for {system}/{date}", fontsize=16, y=1.02)

    # Get data files
    hills_dfs = [] 
    for i in range(n_replicas):
        hills_file = get_file_by_extension(directory, f'{i}.hills')
        hills_dfs.append(read_hills_file(hills_file))
    
    colvar_dfs = []
    for i in range(n_replicas):
        colvar_file = get_file_by_extension(directory, f'{i}.colvar')
        plumed_file = get_file_by_extension(directory, f'plumed_{i}.dat')
        colvar_dfs.append(read_colvar_file(colvar_file, plumed_file))
    
    # Plot hills in first row
    for i in range(n_replicas):
        ax = axs[0, i]
        ax.plot(hills_dfs[i]['time'], hills_dfs[i]['height'])
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel("Height [kJ/mol]")
        ax.set_title(f'Hills {i+1}')

    # Plot first CV trajectory in second row
    for i in range(n_replicas):
        ax = axs[1, i]
        ax.plot(colvar_dfs[i]["time"] * TIMESTEP * STRIDE, colvar_dfs[i][cvs[0]])
        ax.set_title(f'{cvs[0]} {i+1}')
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel(resolve_cv_label(cvs[0]))

    # Plot second CV trajectory in third row
    for i in range(n_replicas):
        ax = axs[2, i]
        ax.plot(colvar_dfs[i]["time"] * TIMESTEP * STRIDE, colvar_dfs[i][cvs[1]])
        ax.set_title(f'{cvs[1]} {i+1}')
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel(resolve_cv_label(cvs[1]))

    # Plot evolution of potential energy and bias in the fourth row
    from src.analysis.traj import plot_biases
    for i in range(n_replicas):
        axs[3, i] = plot_biases(colvar_dfs[i], axs[3, i], timestep=TIMESTEP, stride=STRIDE)

    # Make all axes square and equal
    for ax in axs.flat:
        ax.set_aspect('auto')
        ax.grid(True, linestyle='--', alpha=0.3)

    # Save with high quality settings
    plt.savefig(
        f"{directory}/{system}_replica_trajectories.png",
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    plt.close()


def reconstruct_dcd_from_nc(nc_file, output_dir):
    nc = netCDF4.Dataset(nc_file, 'r')
    positions = nc.variables['positions']

    print(nc.variables['positions'].shape)
    assert 0 == 1

if __name__ == "__main__":
    # directory = "/home/jakub/phd/openmm-md/data/241010_FoldingUponBinding/output/CD28-G/241213-ReplicaPBC"
    directory = "../data/241010_FoldingUponBinding/output/CD28-G/241216-ReplicaPBC"
    system = "CD28-G"
    cvs = ["cmap", "d"]
    n_replicas = 4
    plot_replica_trajectory_analysis(directory, system, cvs, n_replicas)
    plot_replica_analysis(directory, system, cvs, n_replicas)
    nc_file = get_file_by_extension(directory, "replica_exchange.nc")
    reconstruct_dcd_from_nc(nc_file, directory)

