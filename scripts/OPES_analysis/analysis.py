import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import logging

from src.analysis.fes import compute_fes
from src.analysis.plot import plot_2d_fes, plot_1d_fes

from src.analysis.fes import load_fes

logger = logging.getLogger(__name__)

TIMESTEP = 0.002 * 10**-3  # in ns
OPES_PACE = 500
TEMP = 300
STRIDE = 500


def consider_walls(df):
    # compine opes.bias + uwall.bias and create a new column total_bias
    df["total_bias"] = df["opes.bias"] + df["uwall.bias"]
    return df

def plot_all_fes(directory, target, binder, num_runs, labels):
    """Load FES data and create plots
    
    Args:
        directory: Base directory containing run folders
        target: Target protein name 
        binder: Binder name
        num_runs: Number of runs to plot
    """
    # TODO: check, this might be broken now!!!!
    fes_data = []
    for run in range(1, num_runs + 1):
        fes_filepath = f"{directory}/{binder}_{run}/{target}_{binder}_fes.h5py"
        _, fes, cv1_bins, cv2_bins = load_fes(fes_filepath)
        fes_data.append((cv1_bins, cv2_bins, fes))

    plot_all_fes_from_data(fes_data, directory,target, binder, labels=labels)
    outfile = f"{directory}/{binder}_all_fes.png"
    plot_all_fes_from_data(fes_data, outfile, target, binder, ['cmap', 'd'], labels)


from src.analysis.traj import plot_biases

def plot_opes_values(colvar_df, axs):
    ax1, ax2, ax3 = axs
    ax1.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.rct"])
    ax1.set_title("opes.rct")
    ax2.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.zed"])
    ax2.set_title("opes.zed")
    ax3.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.nker"])
    ax3.set_title("opes.nker")
    return ax1, ax2, ax3

from src.analysis.colvar import plot_colvar_trajectories

def plot_summary(directory, system, simulation_type):
    # Create figure with constrained layout for better spacing
    fig = plt.figure(figsize=(12, 10), constrained_layout=True)
    
    # Create GridSpec to ensure equal-sized, square subplots
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1])
    axs = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(3)])
    
    # Add title with padding
    fig.suptitle(f"Analysis for {system}", fontsize=16, y=1.02)

    cvs, fes, cv1_bins, cv2_bins = load_fes(f"{directory}/{system}_fes.h5py")
    colvar_file = get_file_by_extension(directory, '.colvar')  
    colvar_df = read_colvar_file(colvar_file)

    if simulation_type == 'opes':
        axs[0, :] = plot_opes_values(colvar_df, axs[0, :])
    else:
        logger.info(f"Plotting METAD analysis for {system}, hence some empty plots will remain empty")
    
    # Plot colvar trajectories in second row
    axs[1, 0], axs[1, 1] = plot_colvar_trajectories(colvar_df, axs[1, :2], timestep=TIMESTEP, stride=STRIDE)
    axs[1, 2] = plot_biases(colvar_df, axs[1, 2], timestep=TIMESTEP, stride=STRIDE)

    # Plot 2D FES in bottom row
    axs[2, 0] = plot_2d_fes(fes, cv1_bins, cv2_bins, cvs, axs[2, 0], levels=100)
    
    # Plot 1D FES in bottom row
    axs[2, 1], axs[2, 2] = plot_1d_fes(fes, cv1_bins, cv2_bins, cvs, axs[2, 1:])

    # Make all axes square and equal
    for ax in axs.flat:
        ax.set_aspect('auto')
        # Add grid for better readability
        # ax.grid(True, linestyle='--', alpha=0.7)
    
    # Save with high quality settings
    plt.savefig(
        f"{directory}/analysis_summary.png",
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    plt.close()

from matplotlib.animation import FuncAnimation, PillowWriter
from tqdm import tqdm


def plot_colvar_traj_in_fes(directory, system):
    """Create an animated plot of the CV trajectory overlaid on the FES."""
    # Setup figure
    fig, ax = plt.subplots()
    
    # Load data
    cvs, fes, cv1_bins, cv2_bins = load_fes(f"{directory}/{system}_fes.h5py")
    colvar_file = get_file_by_extension(directory, '.colvar')  
    colvar_df = read_colvar_file(colvar_file)
    
    # Get trajectory data (subsample)
    cv1_traj = colvar_df[cvs[0]].values[::1000]
    cv2_traj = colvar_df[cvs[1]].values[::1000]
    
    # Plot FES
    cntr = ax.contourf(cv1_bins, cv2_bins, fes, levels=100, cmap=plt.cm.jet)
    plt.colorbar(cntr, ax=ax, label="FES [kJ/mol]")
    
    # Initialize trajectory line and point
    line, = ax.plot([], [], 'k-', alpha=0.8, linewidth=1)
    point, = ax.plot([], [], 'ko', markersize=8)
    
    ax.set_xlabel(cvs[0])
    ax.set_ylabel(cvs[1])
    ax.set_title(f"System: {system}")
    
    def init():
        line.set_data([], [])
        point.set_data([], [])
        return line, point

    def animate(frame):
        # Update trajectory line
        line.set_data(cv1_traj[:frame], cv2_traj[:frame])
        # Update current point
        point.set_data([cv1_traj[frame]], [cv2_traj[frame]])
        return line, point

    # Create animation with progress bar
    frames = tqdm(range(len(cv1_traj)), desc="Generating animation")
    anim = FuncAnimation(
        fig, 
        animate, 
        init_func=init,
        frames=frames, 
        interval=20,
        blit=True
    )
    
    # Adjust layout and save
    plt.tight_layout()
    writer = PillowWriter(fps=30)
    anim.save(f"{directory}/{system}_trajectory.gif", writer=writer)
    logger.info("Saved trajectory animation")
    plt.close()

def determine_simulation_type(directory):
    plumed_file = get_file_by_extension(directory, 'plumed.dat')
    with open(plumed_file, 'r') as file:
        for line in file:
            if 'OPES' in line:
                logger.info(f"Found OPES in {plumed_file}, assuming OPES simulation")
                return 'opes'
            elif 'METAD' in line:
                logger.info(f"Found METAD in {plumed_file}, assuming METAD simulation")
                return 'metad'
    raise ValueError(f"No simulation type found in {plumed_file}")

from src.analysis.io import read_colvar_file
from src.analysis.kernels import get_sigma
from src.analysis.utils import get_file_by_extension
import pandas as pd

def read_hills_file(hills_file):
    # get first line of kernels file
    with open(hills_file, 'r') as file:
        labels = file.readline().split()[2:]

    df = pd.read_table(
        hills_file,
        dtype=float,
        sep=r"\s+",
        comment="#",
        header=None,
        names=labels
    )
    return df

from src.analysis.utils import fix_fucked_up_naming
def run(project, system, date, recompute, collect_plots):
    date = fix_fucked_up_naming(project, system, date)

    logger.info(f"Processing {project=} {system=} for {date=}")
    directory = f"../../data/{project}/output/{system}/{date}"

    if not os.path.exists(directory):
        logger.error(f"System {system} does not exist for experiment {date}")
        logger.error(f"There is no directory: {directory}")
        return

    colvar_file = get_file_by_extension(directory, '.colvar')  
    colvar_df = read_colvar_file(colvar_file)

    from src.analysis.traj import plot_trajectory
    plot_trajectory(colvar_df, directory, system)

    simulation_type = determine_simulation_type(directory)

    if recompute or not get_file_by_extension(directory, '.h5py', assert_exists=False):
        if simulation_type == 'opes':
            # OPES FES
            cv1_bins, cv2_bins, fes = compute_fes(
                colvar_df,
                sigmas=[get_sigma(directory, cv) for cv in ['cmap', 'd']], 
                temp=300,
                cvs=['cmap', 'd'], 
                outfile=f"{directory}/{system}_fes.h5", 
                bias=['opes.bias', 'uwall.bias'],
                simulation_type='opes'
                )
        elif simulation_type == 'metad':
            # METAD FES
            # cv1_bins, cv2_bins, fes = sum_hills(
            #     colvar_df,
            #     directory=directory,
            #     outfile=f"{directory}/{system}_fes.h5py",
            #     cvs=['cmap', 'd'],
            # )
            cv1_bins, cv2_bins, fes = compute_fes(
                colvar_df,
                sigmas=[get_sigma(directory, cv) for cv in ['cmap', 'd']], 
                temp=300,
                cvs=['cmap', 'd'], 
                outfile=f"{directory}/{system}_fes.h5", 
                bias=['metad.bias', 'uwall.bias'],
                simulation_type='metad'
                )
        else:
            raise ValueError(f"Unknown simulation type: {simulation_type}")

    else:
        cvs, fes, cv1_bins, cv2_bins = load_fes(f"{directory}/{system}_fes.h5py")

    # from src.analysis.plot import plot_2d_fes
    # fig, ax = plt.subplots()
    # plot_2d_fes(fes, cv1_bins, cv2_bins, cvs, ax, levels=100)
    # plt.savefig(f"{directory}/{system}_fes.png", dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    # plt.close()

    plot_colvar_traj_in_fes(directory, system)

    plot_summary(
        directory, 
        system,
        simulation_type
    )

    # if collect_plots:
    #     save_dir = f"../../data/241010_FoldingUponBinding/output/{date}/plots"
    #     os.makedirs(save_dir, exist_ok=True)
    #     shutil.copy(f"{directory}/analysis_summary.png", f"{save_dir}/{target}_{binder}_{run}.png")

    
    # if collect_plots:        
    #     system_directory = "/".join(directory.split("/")[:-1])
            
    #     plot_colvar_traj_in_fes(system_directory, target, binder, num_runs)
    #     shutil.copy(f"{system_directory}/{target}_{binder}_all_trajectories.gif", f"{save_dir}/{target}_{binder}_all_trajectories.gif")
    #     plot_all_fes(system_directory, target, binder, num_runs, labels=[f"{run=}" for run in range(1, num_runs + 1)])
    #     shutil.copy(f"{system_directory}/{binder}_all_fes.png", f"{save_dir}/{target}_{binder}_all_fes.png")

        

        




def main(project, system, date, recompute=False):
    collect_plots = True
    run(project, system, date, recompute, collect_plots)
    
import fire
if __name__ == "__main__":
    fire.Fire(main)
