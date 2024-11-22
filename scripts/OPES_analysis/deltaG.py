from src.analysis.utils import get_file_by_extension
from src.analysis.colvar import read_colvar_file
from src.analysis.deltaG import compute_deltaG, marginalize_fes
from src.analysis.fes import load_fes
from src.constants import kB, mol
from src.analysis.utils import convert_deltaG_to_kBT
from src.analysis.kernels import get_sigma
from src.analysis.fes import compute_fes
import tempfile
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

TIMESTEP = 0.002
STRIDE = 500
TEMP = 300
SPLIT_RATIO = 10


def compute_deltaG_from_trajectory(directory, colvar_df, threshold):
    split_size = len(colvar_df) // SPLIT_RATIO

    deltaGs = []
    with tempfile.TemporaryDirectory() as temp_dir:
        for i in tqdm(range(1, SPLIT_RATIO + 1), desc=f"Computing FES for chunks"):
            cv1_bins, cv2_bins, fes = compute_fes(
                colvar_df[:int(i*split_size)],
                sigmas=[get_sigma(directory, cv) for cv in ['cmap', 'd']], 
                temp=300,
                cvs=['cmap', 'd'], 
                outfile=f"{temp_dir}/{i}_fes.h5py", 
                bias=['opes.bias', 'uwall.bias']
            )
            fes_1d = marginalize_fes(fes, kBT=kB*TEMP, axis=1)
            deltaG = compute_deltaG(fes_1d, cv2_bins, kB*TEMP, threshold)
            deltaGs.append(deltaG)
    plot_deltaGs(deltaGs, split_size, directory)
    return deltaGs

def plot_deltaGs(deltaGs, split_size, directory):
    # Calculate time points and statistics once
    times = np.arange(len(deltaGs)) * split_size * TIMESTEP * STRIDE
    mean_deltaG = np.mean(deltaGs)
    std_deltaG = np.std(deltaGs)
    
    # Create figure with twin axes
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    
    # Plot on first axis (kJ/mol)
    ax1.scatter(times, deltaGs, marker='x', color='blue', label='ΔG values')
    ax1.axhline(mean_deltaG, color="green", linestyle="--", 
                label=f'Average: {mean_deltaG:.2f} ± {std_deltaG:.2f} kJ/mol')
    
    # Plot on second axis (kBT)
    deltaGs_kbt = [convert_deltaG_to_kBT(dg, TEMP) for dg in deltaGs]
    ax2.scatter(times, deltaGs_kbt, alpha=0)  # invisible points to set scale
    
    # Labels and styling
    ax1.set_xlabel('Time [ps]')
    ax1.set_ylabel('ΔG [kJ/mol]')
    ax2.set_ylabel('ΔG [kBT]')
    
    plt.title('Binding Affinity Estimate from Cumulative Trajectory')
    
    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(lines1, labels1, loc='upper right')
    
    # Save and close
    plt.tight_layout()
    plt.savefig(f"{directory}/deltaGs.png", dpi=300, bbox_inches='tight')
    plt.close()

from src.analysis.utils import fix_fucked_up_naming
def run(project, system, date):
    date = fix_fucked_up_naming(project, system, date)
    directory = f"../../data/{project}/output/{system}/{date}"
    # 1. Get the binding affinity from the FES
    cvs, fes, cv1_bins, cv2_bins = load_fes(f"{directory}/{system}_fes.h5py")
    assert cvs[1] == 'd', "Second CV is not distance between COMs, we do not support other CVs yet"
    fes_1d = marginalize_fes(fes, kBT=kB*TEMP, axis=1)
    deltaG = compute_deltaG(fes_1d, cv2_bins, kB*TEMP, threshold=1.5)
    print(f"Delta G [kJ/mol]: {deltaG:.2f} kJ/mol")
    print(f"Delta G [kBT]: {convert_deltaG_to_kBT(deltaG, TEMP):.2f} kBT")
    
    # 2. Compute the binding affinity for subsets of the trajectory
    colvar_file = get_file_by_extension(directory, ".colvar")
    colvar_df = read_colvar_file(colvar_file)
    deltaGs = compute_deltaG_from_trajectory(directory, colvar_df[:int(0.1*len(colvar_df))], threshold=1.5)
    


def main(project, system, date):
    run(project, system, date)
    
import fire
if __name__ == "__main__":
    fire.Fire(main)
