import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import h5py
import re
import logging

from src.analysis.fes import compute_fes
from src.constants import kB
import pandas as pd
from src.analysis.deltaG import marginalize_fes
    

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

def get_sigma(directory, target, binder, run, cv):
    # get the sigma from the kernels file
    kernels_file = f"{directory}/{target}_{binder}.kernels"
    df = pd.read_table(
        kernels_file,
        dtype=float,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["time", "cmap", "d", "sigma_cmap", "sigma_d", "height", "logweight"]
    )
    sigma = df[f"sigma_{cv}"].iloc[-1]
    return sigma.item()


from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_2d_fes(fes, cv1_bins, cv2_bins, cvs, ax):
    cv1, cv2 = cvs

    ax.set_xlabel(cv1)
    ax.set_ylabel(cv2)

    # Add filled contours
    cntr = ax.contourf(cv1_bins, cv2_bins, fes, levels=range(0, 120, 5), cmap=plt.cm.jet)

    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad="20%")
    cbar = plt.colorbar(cntr, cax=cax, label="FES [kJ/mol]")
    cbar.ax.yaxis.set_label_position("left")

    # Add kBT scale
    cbar_kBT = cbar.ax.twinx()
    cbar_kBT.set_ylim(cbar.vmin / (kB * TEMP), cbar.vmax / (kB * TEMP))
    cbar_kBT.set_ylabel("FES [kBT]", rotation=270, labelpad=15)

    return ax

def get_sigmas(directory, target, binder, run, cvs):
    assert len(cvs) == 2, "Only 2D FES are supported"
    return [get_sigma(directory, target, binder, run, cv) for cv in cvs]

# compute FES
def read_colvar_file(filename):
    # Read first line to get column names
    with open(filename) as f:
        header_line = f.readline().strip()
    
    # Parse column names from FIELDS line
    if not header_line.startswith('#! FIELDS'):
        raise ValueError("First line must start with '#! FIELDS'")
    column_names = header_line.replace('#! FIELDS', '').strip().split()
    
    # Read data using column names from file
    colvar_df = pd.read_table(
        filename,
        dtype=float,
        sep=r"\s+",
        comment="#", 
        header=None,
        names=column_names
    )

    # remove the unbiased bit
    # PLUMED:   adaptive SIGMA will be used, with ADAPTIVE_SIGMA_STRIDE = 5000
    # PLUMED:     thus the first x kernel depositions will be skipped, x = ADAPTIVE_SIGMA_STRIDE/PACE = 10
    # get the pace from the _plumed.dat file
    filename = filename.replace(".colvar", "_plumed.dat")
    with open(filename, "r") as f:
        content = f.read()
    
    # depositing kernels
    pace_pattern = r"opes:\s*OPES_METAD[\s\S]*?PACE=(\d+)"
    match = re.search(pace_pattern, content)
    if match:
        pace = int(match.group(1))
        logger.info(f"Pace value: {pace}")
    else:
        raise ValueError("PACE not found in the file")
    
    # printing to COLVAR file
    stride_pattern = r"PRINT.*STRIDE=(\d+)"
    match = re.search(stride_pattern, content)
    if match:
        stride = int(match.group(1))
        logger.info(f"Stride value: {stride}")
    else:
        raise ValueError("Stride not found in the file")
    
    # assume first x = 10 kernel depositions are the ones to skip
    # compute how many printing steps correspond to that
    num_steps = pace * 10
    strides_to_skip = num_steps // stride
    logger.info(f"Skipping {strides_to_skip} steps")

    # look at the opes.bias in the strides_to_skip steps and assert they are constant
    assert colvar_df["opes.bias"].iloc[:strides_to_skip].std() < 1e-6, "opes.bias is not constant in non-biased part"
    colvar_df = colvar_df.iloc[strides_to_skip:]
    return colvar_df

def plot_1d_fes(fes, cv1_bins, cv2_bins, cvs, axs):
    ax1, ax2 = axs
    cv1, cv2 = cvs
    # plot along first CV
    # axis = 0, as we keep the first CV
    fes_1d = marginalize_fes(fes, kBT=kB*TEMP, axis=0)
    ax1.plot(cv1_bins, fes_1d)
    ax1.set_xlabel(cv1)
    ax1.set_ylabel("FES [kJ/mol]")

    # plot along second CV
    fes_1d = marginalize_fes(fes, kBT=kB*TEMP, axis=1)
    ax2.plot(cv2_bins, fes_1d)
    ax2.set_xlabel(cv2)
    ax2.set_ylabel("FES [kJ/mol]")
    ax2.invert_xaxis()
    return ax1, ax2

def load_fes(filepath: str):
    with h5py.File(filepath, "r") as f:
        cvs = f.attrs['cvs']
        fes = f["fes"][:]
        cv1_bins = f[f"{cvs[0]}_bins"][:]
        cv2_bins = f[f"{cvs[1]}_bins"][:]
        first_axis = f.attrs['1st_axis']
        second_axis = f.attrs['2nd_axis']

    assert first_axis in cvs and second_axis in cvs, "First and second axis must be in CVs"
    assert first_axis != second_axis, "First and second axis must be different"
    assert first_axis == cvs[0] and second_axis == cvs[1], "First axis must be the first CV and second axis must be the second CV"
    assert fes.shape[0] == len(cv1_bins) and fes.shape[1] == len(cv2_bins), "FES must have the same shape as the CV bins"
    return cvs, fes, cv1_bins, cv2_bins

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


def plot_all_fes(directory, target, binder, num_runs):
    # plot all the FES for a given system in a single plot
    # make it scalable so that you have the same number of columns as runs
    # first row: full FES
    # second row: 1D FES along first CV
    # third row: 1D FES along second CV

    # Calculate figure width based on number of runs to maintain square subplots
    # Height is fixed at 15 inches, width scales with number of runs
    fig_height = 15
    subplot_size = fig_height / 3  # Each subplot should be square
    fig_width = subplot_size * num_runs * 1.2  # Add 20% more width for spacing
    
    # Create figure with square-looking subplots
    fig = plt.figure(figsize=(fig_width, fig_height))
    fig.suptitle(f"{binder=} and {target=}", fontsize=16, y=0.98)

    # Use gridspec with adjusted spacing
    gs = fig.add_gridspec(3, num_runs, wspace=0.3, hspace=0.3, top=0.95)
    axs = np.empty((3, num_runs), dtype=object)
    
    # Create subplots using gridspec
    for i in range(3):
        for j in range(num_runs):
            axs[i,j] = fig.add_subplot(gs[i,j])
            if i == 0:  # Add run number at the top of each column
                axs[i,j].set_title(f"run={j+1}", pad=5)  # Reduced pad value

    cv1_mins = []
    cv1_maxs = []
    cv2_mins = []
    cv2_maxs = []

    for run in range(1, num_runs + 1):
        fes_filepath = f"{directory}/{binder}_{run}/{target}_{binder}_fes.h5py"
        cvs, fes, cv1_bins, cv2_bins = load_fes(fes_filepath)
        axs[0, run - 1] = plot_2d_fes(fes, cv1_bins, cv2_bins, cvs, axs[0, run - 1])
        axs[1, run - 1], axs[2, run - 1] = plot_1d_fes(fes, cv1_bins, cv2_bins, cvs, axs[1:3, run - 1])

        cv1_mins.append(cv1_bins[0])
        cv1_maxs.append(cv1_bins[-1])
        cv2_mins.append(cv2_bins[0])
        cv2_maxs.append(cv2_bins[-1])

    cv1_min = min(cv1_mins)
    cv1_max = max(cv1_maxs)
    cv2_min = min(cv2_mins)
    cv2_max = max(cv2_maxs)

    # set limits on all subplots
    # first row, use the cv1 and so for each of the axis
    for j in range(num_runs):
        axs[0, j].set_xlim(cv1_min, cv1_max)
        axs[0, j].set_ylim(cv2_min, cv2_max)
    # second row, change based on cv1
    for j in range(num_runs):
        axs[1, j].set_xlim(cv1_min, cv1_max)
    # third row, change based on cv2
    for j in range(num_runs):
        # this one is flipped
        axs[2, j].set_xlim(cv2_max, cv2_min)
    


    plt.savefig(
        f"{directory}/{binder}_all_fes.png",
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    plt.close()

    return    

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

            plot_all_fes(system_directory, target, binder, num_runs)
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
