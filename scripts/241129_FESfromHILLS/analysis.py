from src.analysis.fes import compute_fes_from_hills
import fire
import os
import numpy as np
import matplotlib.pyplot as plt
from src.analysis.plot import plot_2d_fes, plot_1d_fes
from src.analysis.utils import get_file_by_extension
from src.analysis.io import read_hills_file
from src.analysis.fes import load_fes

def main(system: str, project: str, date: str):

    directory = f"../../data/{project}/output/{system}/{date}"
    hills_file = get_file_by_extension(
        directory=directory,
        extension=".hills"
    )

    hills_df = read_hills_file(hills_file)

    # assert hills_df['biasf'] is all the same value
    assert hills_df['biasf'].std() < 1e-6, "biasf is not constant"

    cvs = ["cmap", "d"]
    if not os.path.exists(f"{directory}/{system}_fes_hills.h5"):
        cv1_bins_hills, cv2_bins_hills, fes_hills = compute_fes_from_hills(
            hills_df=hills_df,
            temp=300,
            cvs=cvs,
            biasfactor=hills_df['biasf'].iloc[0],
            outfile=f"{directory}/{system}_fes_hills.h5",
            n_bins=100
        )
    else:
        _, fes_hills, cv1_bins_hills, cv2_bins_hills = load_fes(f"{directory}/{system}_fes_hills.h5")

    try:
        fes_reweighted_file = get_file_by_extension(
            directory=directory,
            extension="reweighted.h5"
        )
        _, fes_reweighted, cv1_bins_reweighted, cv2_bins_reweighted = load_fes(fes_reweighted_file)
    except FileNotFoundError:
        from src.analysis.io import read_colvar_file
        from src.analysis.kernels import get_sigma
        from src.analysis.fes import compute_fes
        colvar_df = read_colvar_file(get_file_by_extension(directory, ".colvar"))
        cv1_bins_reweighted, cv2_bins_reweighted, fes_reweighted = compute_fes(
            colvar_df=colvar_df,
            sigmas=[get_sigma(directory, cv) for cv in cvs],
            temp=300,
            cvs=cvs,
            outfile=f"{directory}/{system}_fes_reweighted.h5",
            n_bins=100,
            simulation_type="metad"
        )
    fes_reweighted = fes_reweighted.T

    # Create figure with constrained layout for better spacing
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    
    # Create GridSpec to ensure equal-sized subplots with extra left margin
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1], left=0.1)
    axs = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(3)])

    # Add main title with padding
    fig.suptitle(f"FES Analysis for {system}", fontsize=16, y=1.02)

    levels = np.unique(np.linspace(fes_hills.min(), fes_hills.max(), 100).astype(np.int32))

    # assert increasing levels
    assert np.all(np.diff(levels) > 0), f"Levels are not increasing, they are {levels}"


    # Plot FES comparisons
    axs[0, 0] = plot_2d_fes(fes_hills, cv1_bins_hills, cv2_bins_hills, cvs, axs[0, 0], levels=levels)
    axs[1, 0] = plot_2d_fes(fes_reweighted, cv1_bins_reweighted, cv2_bins_reweighted, cvs, axs[1, 0], levels=levels)
    axs[0, 1], axs[0, 2] = plot_1d_fes(fes_hills, cv1_bins_hills, cv2_bins_hills, cvs, axs[0, 1:])
    axs[1, 1], axs[1, 2] = plot_1d_fes(fes_reweighted, cv1_bins_reweighted, cv2_bins_reweighted, cvs, axs[1, 1:])

    # Sync y-axis limits between corresponding 1D plots
    axs[1, 1].set_ylim(axs[0, 1].get_ylim())
    axs[1, 2].set_ylim(axs[0, 2].get_ylim())

    axs[2, 2].set_ylim(axs[0, 2].get_ylim())
    axs[2, 1].set_ylim(axs[0, 1].get_ylim())

    # in the third row, plot the 1d fes from both on the same plot
    # get the values from the first and second row
    hills_line = axs[0, 1].get_lines()[0]  # get first line from the plot
    reweighted_line = axs[1, 1].get_lines()[0]
    
    x_hills, y_hills = hills_line.get_data()
    x_reweighted, y_reweighted = reweighted_line.get_data()
    axs[2, 1].plot(x_hills, y_hills, label='Summed Hills')
    axs[2, 1].plot(x_reweighted, y_reweighted, label='Reweighting')
    axs[2, 1].legend()

    # do the same for the third column
    hills_line = axs[0, 2].get_lines()[0]  # get first line from the plot
    reweighted_line = axs[1, 2].get_lines()[0]
    
    x_hills, y_hills = hills_line.get_data()
    x_reweighted, y_reweighted = reweighted_line.get_data()
    axs[2, 2].plot(x_hills, y_hills, label='Summed Hills')
    axs[2, 2].plot(x_reweighted, y_reweighted, label='Reweighting')
    axs[2, 2].legend()

    # Add row labels with more space on left
    fig.text(0.05, 0.75, 'Summed Hills', ha='left', va='center', fontsize=12, rotation=90)
    fig.text(0.05, 0.25, 'Reweighting', ha='left', va='center', fontsize=12, rotation=90)

    # Make all axes square and equal
    for ax in axs.flat:
        ax.set_aspect('auto')

    # Save with high quality settings
    plt.savefig(
        f"{directory}/{system}_fes_comparison.png",
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    plt.close()

if __name__ == "__main__":
    fire.Fire(main)
