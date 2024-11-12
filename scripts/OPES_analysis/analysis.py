import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import h5py
import re
import logging

from src.analysis.fes import compute_fes
from src.constants import kB
from src.analysis.colvar import read_colvar_file
from src.analysis.kernels import get_sigmas
    

logger = logging.getLogger(__name__)

TIMESTEP = 0.002 * 10**-3  # in ns
OPES_PACE = 500
TEMP = 300


def plot_colvar_trajectories(df, axs):
    ax1, ax2 = axs

    # skip every 100th point, to smooth out the trajectories
    df = df.iloc[::100]

    # TODO: make this more general, cmap and d are hard-locked
    ax1.plot(df["time"] * TIMESTEP * OPES_PACE, df["cmap"], label='CMAP')
    ax2.plot(df["time"] * TIMESTEP * OPES_PACE, df["d"], label='Distance')

    ax1.set_title("CMAP")
    ax1.set_xlabel("Time [ns]")
    ax1.set_ylabel("CMAP [#]")
    ax1.legend()

    ax2.set_title("Distance")
    ax2.set_xlabel("Time [ns]")
    ax2.set_ylabel("Distance [nm]")
    ax2.legend()
    return ax1, ax2


def consider_walls(df):
    # compine opes.bias + uwall.bias and create a new column total_bias
    df["total_bias"] = df["opes.bias"] + df["uwall.bias"]
    return df

from src.analysis.plot import plot_all_fes_from_data
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

from src.analysis.plot import plot_2d_fes, plot_1d_fes
from src.analysis.fes import load_fes
# compute FES

def plot_bias(colvar_df, ax):
    ax.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.bias"] + colvar_df["uwall.bias"], label='total')
    ax.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.bias"], label='opes')
    ax.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["uwall.bias"], label='uwall')
    ax.legend()
    ax.set_title("Bias")
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Bias [kJ/mol]")
    return ax

def plot_opes_values(colvar_df, axs):
    ax1, ax2, ax3 = axs
    ax1.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.rct"])
    ax1.set_title("opes.rct")
    ax2.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.zed"])
    ax2.set_title("opes.zed")
    ax3.plot(colvar_df["time"] * TIMESTEP * OPES_PACE, colvar_df["opes.nker"])
    ax3.set_title("opes.nker")
    return ax1, ax2, ax3

def plot_everything(directory, target, binder, run):
    # Create figure with constrained layout for better spacing
    fig = plt.figure(figsize=(12, 10), constrained_layout=True)
    
    # Create GridSpec to ensure equal-sized, square subplots
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1])
    axs = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(3)])
    
    # Add title with padding
    fig.suptitle(f"Analysis for {target}_{binder} run {run}", fontsize=16, y=1.02)

    cvs, fes, cv1_bins, cv2_bins = load_fes(f"{directory}/{target}_{binder}_fes.h5py")
    colvar_df = read_colvar_file(f"{directory}/{target}_{binder}.colvar")

    axs[0, :] = plot_opes_values(colvar_df, axs[0, :])
    
    # Plot colvar trajectories in second row
    axs[1, 0], axs[1, 1] = plot_colvar_trajectories(colvar_df, axs[1, :2])
    axs[1, 2] = plot_bias(colvar_df, axs[1, 2])

    # Plot 2D FES in bottom row
    axs[2, 0] = plot_2d_fes(fes, cv1_bins, cv2_bins, cvs, axs[2, 0])
    
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

from concurrent.futures import ThreadPoolExecutor
from functools import partial

def process_single_run(run, directory, target, binder, shared_fig, shared_axes):
    cvs, fes, cv1_bins, cv2_bins = load_fes(f"{directory}/{binder}_{run}/{target}_{binder}_fes.h5py")
    colvar_df = read_colvar_file(f"{directory}/{binder}_{run}/{target}_{binder}.colvar")
    
    # Get trajectory data (subsample)
    cv1_traj = colvar_df[cvs[0]].values[::1000]
    cv2_traj = colvar_df[cvs[1]].values[::1000]
    
    ax = shared_axes[run-1]
    
    # Plot FES
    cntr = ax.contourf(cv1_bins, cv2_bins, fes, levels=range(0, 120, 5), cmap=plt.cm.jet)
    plt.colorbar(cntr, ax=ax, label="FES [kJ/mol]")
    
    # Initialize trajectory line and point
    line, = ax.plot([], [], 'k-', alpha=0.8, linewidth=1)
    point, = ax.plot([], [], 'ko', markersize=8)
    
    ax.set_xlabel(cvs[0])
    ax.set_ylabel(cvs[1])
    ax.set_title(f"Run {run}")
    
    return {
        'ax': ax,
        'line': line,
        'point': point,
        'cv1_traj': cv1_traj,
        'cv2_traj': cv2_traj
    }

def plot_colvar_traj_in_fes(directory, target, binder, num_runs):
    # Create a wide figure to accommodate all runs
    fig_width = 6 * num_runs  # 6 inches per subplot
    fig, axes = plt.subplots(1, num_runs, figsize=(fig_width, 6))
    if num_runs == 1:
        axes = [axes]
    
    # Process all runs in parallel
    with ThreadPoolExecutor() as executor:
        process_run = partial(process_single_run, 
                            directory=directory, 
                            target=target, 
                            binder=binder, 
                            shared_fig=fig, 
                            shared_axes=axes)
        
        run_data = list(executor.map(process_run, range(1, num_runs + 1)))
    
    # Find the maximum trajectory length
    max_frames = max(len(data['cv1_traj']) for data in run_data)
    
    def init():
        elements = []
        for data in run_data:
            data['line'].set_data([], [])
            data['point'].set_data([], [])
            elements.extend([data['line'], data['point']])
        return elements

    def animate(frame):
        elements = []
        for data in run_data:
            # Handle different trajectory lengths
            curr_frame = min(frame, len(data['cv1_traj'])-1)
            
            # Update trajectory line
            data['line'].set_data(data['cv1_traj'][:curr_frame], 
                                data['cv2_traj'][:curr_frame])
            # Update current point
            data['point'].set_data([data['cv1_traj'][curr_frame]], 
                                 [data['cv2_traj'][curr_frame]])
            elements.extend([data['line'], data['point']])
        return elements

    # Create animation with progress bar
    frames = tqdm(range(max_frames), desc="Generating animation")
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
    anim.save(f"{directory}/{target}_{binder}_all_trajectories.gif", writer=writer)
    logger.info("Saved combined trajectory animation")
    plt.close()


def run(date, systems, num_runs, recompute, collect_plots):
    for system in systems:
        target, binder = system.split("_")
        for run in range(1, num_runs + 1):
            logger.info(f"Processing {system} run {run}")
            directory = f"../../data/241010_FoldingUponBinding/output/{date}/{target}/{binder}_{run}"
            # check system exists, if not write it
            if not os.path.exists(directory):
                logger.warning(f"System {system} does not exist for experiment {date}")
                continue

            if target == "SUMO":
                target = "sumo"
            
            colvar_df = read_colvar_file(f"{directory}/{target}_{binder}.colvar")
            if recompute or not os.path.exists(f"{directory}/{target}_{binder}_fes.h5py"):
                _, _, _ = compute_fes(
                    colvar_df,
                    sigma=get_sigmas(directory, target, binder, run, ['cmap', 'd']), 
                    temp=300, 
                    cvs=['cmap', 'd'], 
                    outfile=f"{directory}/{target}_{binder}_fes.h5py", 
                    bias=['opes.bias', 'uwall.bias']
                )

            plot_everything(directory, target, binder, run)

            if collect_plots:
                save_dir = f"../../data/241010_FoldingUponBinding/output/{date}/plots"
                os.makedirs(save_dir, exist_ok=True)
                shutil.copy(f"{directory}/analysis_summary.png", f"{save_dir}/{target}_{binder}_{run}.png")

        
        if collect_plots:        
            system_directory = "/".join(directory.split("/")[:-1])
                
            plot_colvar_traj_in_fes(system_directory, target, binder, num_runs)
            shutil.copy(f"{system_directory}/{target}_{binder}_all_trajectories.gif", f"{save_dir}/{target}_{binder}_all_trajectories.gif")
            plot_all_fes(system_directory, target, binder, num_runs, labels=[f"{run=}" for run in range(1, num_runs + 1)])
            shutil.copy(f"{system_directory}/{binder}_all_fes.png", f"{save_dir}/{target}_{binder}_all_fes.png")

        




def main():
    global num_runs
    systems = ['A-synuclein_alpha', 'A-synuclein_general', 'CD28_alpha', 'CD28_beta', 'CD28_partial', 'CD28_general']
    date = "241029"
    num_runs = 5
    recompute = False
    collect_plots = True
    run(date, systems, num_runs, recompute, collect_plots)
    date = "241028"
    systems = ['p53_1', 'p53_2', 'p53_end']
    run(date, systems, num_runs, recompute, collect_plots)
    # systems = ['SUMO_1', 'SUMO_1c']
    # run(date, systems, num_runs, recompute, collect_plots)
   

if __name__ == "__main__":
    main()
