import os
import plumed
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import h5py

import logging

from src.analysis.fes import compute_fes
import pandas as pd

logger = logging.getLogger(__name__)

TIMESTEP = 0.002 * 10**-3  # in ns
OPES_PACE = 500


def plot_colvar_trajectories(directory, target, binder, run):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Colvar Trajectories for {target=} {binder=} {run=}")

    # get the colvar trajectory
    colvar_file = f"{directory}/{target}_{binder}.colvar"
    try:
        df = read_colvar_file(colvar_file)
    except FileNotFoundError:
        # If the file doesn't exist, continue to the next iteration
        print(f"Colvar file {colvar_file} does not exist")
        return

    ax1.plot(df["time"] * TIMESTEP * OPES_PACE, df["cmap"], label=binder)
    ax2.plot(df["time"] * TIMESTEP * OPES_PACE, df["d"], label=binder)

    ax1.set_title("CMAP")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("CMAP")
    ax1.legend()

    ax2.set_title("Distance")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Distance")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(f"{directory}/colvar_trajectories.png", dpi=300)


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


# def compute_FES(directory, target, binder, run, recompute=False):
#     colvar_file = f"{directory}/{target}_{binder}.colvar"
#     df = plumed.read_as_pandas(colvar_file)
#     df = consider_walls(df)
#     # save df into file as colvar
#     plumed.write_pandas(df, f"{directory}/{target}_{binder}_with_walls.colvar")

#     sigma_cmap = get_sigma(directory, target, binder, run, "cmap")
#     sigma_d = get_sigma(directory, target, binder, run, "d")

#     # if the output file exists, skip the computation
#     if not os.path.exists(f"{directory}/{target}_{binder}_fes.dat") or recompute:
#         if recompute:
#             logger.info(f"Recomputing FES for {target} {binder} {run}")
#         else:
#             logger.info(f"Computing FES for {target} {binder} {run}")
#         # run the script FES_from_Reweighting.py
#         subprocess.run(
#             [
#                 "python",
#                 "FES_from_Reweighting.py",
#                 "--colvar",
#                 f"{directory}/{target}_{binder}_with_walls.colvar",
#                 "--outfile",
#                 f"{directory}/{target}_{binder}_fes.dat",
#                 "--sigma",
#                 f"{sigma_cmap},{sigma_d}",
#                 "--temp",
#                 "300",
#                 "--bias",
#                 "opes.bias",
#                 "--cv",
#                 "cmap,d",
#             ]
#         )

#     # plot the 2d FES
#     plot_2d_fes(directory, target, binder, run)


from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_2d_fes(directory, target, binder, run):
    # Read FES data
    with h5py.File(f"{directory}/{target}_{binder}_fes.h5py", "r") as f:
        fes = f["fes"][:]
        cmap = f["cmap_bins"][:]
        d = f["d_bins"][:]

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))  # Add figure size for better proportions

    # Add contour lines
    ax.contour(cmap, d, fes, levels=range(0, 120, 5), linewidths=0.5, colors="k")
    ax.set_xlabel("Contact Map [# of contacts]")
    ax.set_ylabel("Distance between COMs [nm]")

    # Add filled contours
    cntr = ax.contourf(cmap, d, fes, levels=range(0, 120, 5), cmap=plt.cm.jet)

    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad="20%")
    cbar = plt.colorbar(cntr, cax=cax, label="FES [kJ/mol]")
    cbar.ax.yaxis.set_label_position("left")

    # Add kBT scale
    from src.constants import kB
    TEMP = 300
    cbar_kBT = cbar.ax.twinx()
    cbar_kBT.set_ylim(cbar.vmin / (kB * TEMP), cbar.vmax / (kB * TEMP))
    cbar_kBT.set_ylabel("FES [kBT]", rotation=270, labelpad=15)

    plt.title(f"{target} {binder} (run {run})")
    plt.savefig(f"{directory}/fes_2d.png", dpi=300)
    plt.close()

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
    return colvar_df

def main():
    # List of all systems
    systems = ['A-synuclein_alpha', 'A-synuclein_general', 'CD28_alpha', 'CD28_beta', 'CD28_partial']
    date = "241029"
    num_runs = 5
    recompute = True

    for system in systems:
        target, binder = system.split("_")
        for run in range(1, num_runs + 1):
            logger.info(f"Processing {system} run {run}")
            directory = f"../../data/241010_FoldingUponBinding/output/{date}/{target}/{binder}_{run}"
            # check system exists, if not write it
            if not os.path.exists(directory):
                logger.warning(f"System {system} does not exist for experiment {date}")
                continue
            plot_colvar_trajectories(directory, target, binder, run)
                
            colvar_df = read_colvar_file(f"{directory}/{target}_{binder}.colvar")
            _, _, _ = compute_fes(
                colvar_df,
                sigma=get_sigmas(directory, target, binder, run, ['cmap', 'd']), 
                temp=300, 
                cvs=['cmap', 'd'], 
                outfile=f"{directory}/{target}_{binder}_fes.h5py", 
                bias=['opes.bias', 'uwall.bias']
            )

            plot_2d_fes(directory, target, binder, run)


if __name__ == "__main__":
    main()
