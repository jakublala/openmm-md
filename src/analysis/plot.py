import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

# HACK:
TIMESTEP = 0.002 * 10**-3  # in ns
OPES_PACE = 500
TEMP = 300
from src.constants import kB
from src.analysis.deltaG import marginalize_fes
from src.analysis.fes import load_fes


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


def plot_all_fes_from_data(fes_data, outfile, target, binder, cvs, labels):
    """Plot all FES from preloaded data
    
    Args:
        fes_data: List of tuples (cv1_bins, cv2_bins, fes) for each run
        target: Target protein name
        binder: Binder name
    """
    num_runs = len(fes_data)
    
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
                axs[i,j].set_title(labels[j], pad=5)

    cv1_mins = []
    cv1_maxs = []
    cv2_mins = []
    cv2_maxs = []

    for run in range(num_runs):
        cv1_bins, cv2_bins, fes = fes_data[run]
        
        axs[0, run] = plot_2d_fes(fes, cv1_bins, cv2_bins, cvs, axs[0, run])
        axs[1, run], axs[2, run] = plot_1d_fes(fes, cv1_bins, cv2_bins, cvs, axs[1:3, run])

        cv1_mins.append(cv1_bins[0])
        cv1_maxs.append(cv1_bins[-1])
        cv2_mins.append(cv2_bins[0])
        cv2_maxs.append(cv2_bins[-1])

    cv1_min = min(cv1_mins)
    cv1_max = max(cv1_maxs)
    cv2_min = min(cv2_mins)
    cv2_max = max(cv2_maxs)

    # set limits on all subplots
    for j in range(num_runs):
        axs[0, j].set_xlim(cv1_min, cv1_max)
        axs[0, j].set_ylim(cv2_min, cv2_max)
        axs[1, j].set_xlim(cv1_min, cv1_max)
        axs[2, j].set_xlim(cv2_max, cv2_min)  # Flipped for third row

    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    plt.close()
    return True
