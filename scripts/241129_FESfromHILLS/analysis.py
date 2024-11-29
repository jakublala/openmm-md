from src.analysis.fes import compute_fes_from_hills
import fire
import os
import numpy as np
import matplotlib.pyplot as plt
from src.analysis.plot import plot_2d_fes, plot_1d_fes
from src.analysis.utils import get_file_by_extension
from src.analysis.io import read_hills_file
from src.analysis.fes import load_fes
def main(system: str = "ASYN-A"):
    directory = f"../../data/241010_FoldingUponBinding/output/{system}/241128-MetaD"
    hills_file = get_file_by_extension(
        directory=directory,
        extension=".hills"
    )

    hills_df = read_hills_file(hills_file)

    # assert hills_df['biasf'] is all the same value
    assert hills_df['biasf'].std() < 1e-6, "biasf is not constant"

    cvs = ["cmap", "d"]
    if not os.path.exists(f"{system}_fes_hills.h5"):
        cv1_bins_hills, cv2_bins_hills, fes_hills = compute_fes_from_hills(
            hills_df=hills_df,
            temp=300,
            cvs=cvs,
            biasfactor=hills_df['biasf'].iloc[0],
            outfile=f"{system}_fes_hills.h5",
            n_bins=100
        )
    else:
        _, fes_hills, cv1_bins_hills, cv2_bins_hills = load_fes(f"{system}_fes_hills.h5")

    fes_reweighted_file = get_file_by_extension(
        directory=directory,
        extension=".h5"
    )
    if not os.path.exists(fes_reweighted_file):
        from src.analysis.io import read_colvar_file
        from src.analysis.kernels import get_sigmas
        colvar_df = read_colvar_file(get_file_by_extension(directory, ".colvar"))
        cv1_bins_reweighted, cv2_bins_reweighted, fes_reweighted = compute_fes(
            colvar_df=colvar_df,
            sigmas=[get_sigma(directory, cv) for cv in cvs],
            temp=300,
            cvs=cvs,
            outfile=fes_reweighted_file,
        )
    _, fes_reweighted, cv1_bins_reweighted, cv2_bins_reweighted = load_fes(fes_reweighted_file)

    fig, axs = plt.subplots(2, 3, figsize=(12, 8))


    levels = np.linspace(fes_hills.min(), fes_hills.max(), 100).astype(np.int32)

    axs[0, 0] = plot_2d_fes(fes_hills, cv1_bins_hills, cv2_bins_hills, cvs, axs[0, 0], levels=levels)
    axs[1, 0] = plot_2d_fes(fes_reweighted, cv1_bins_reweighted, cv2_bins_reweighted, cvs, axs[1, 0], levels=levels)
    axs[0, 1], axs[0, 2] = plot_1d_fes(fes_hills, cv1_bins_hills, cv2_bins_hills, cvs, axs[0, 1:])
    axs[1, 1], axs[1, 2] = plot_1d_fes(fes_reweighted, cv1_bins_reweighted, cv2_bins_reweighted, cvs, axs[1, 1:])

    # Add rotated titles for each row
    fig.text(0.4, 0.95, 'Summed Hills', ha='center', fontsize=12)
    fig.text(0.4, 0.45, 'Reweighting', ha='center', fontsize=12)


    # fix the same xticks, yticks for both same plots
    plt.savefig(f"{system}_fes_comparison.png", dpi=300)

if __name__ == "__main__":
    fire.Fire(main)
