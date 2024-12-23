import matplotlib.pyplot as plt
from src.analysis.plot import plot_1d_fes
import h5py
from src.analysis.fes import load_fes

def main():
    fes_filepath = "fes.h5"

    # load the file
    cvs, fes, cv1_bins, cv2_bins = load_fes(fes_filepath)

    from src.constants import kB
    from src.analysis.deltaG import marginalize_fes
    TEMP = 300

    # Create figure and axes properly
    fig, ax = plt.subplots()  # This is the correct way to create figure and axes
    
    fes_1d = marginalize_fes(fes, kBT=kB*TEMP, axis=0, dx=cv2_bins[1]-cv2_bins[0])
    ax.plot(cv2_bins, fes_1d)
    ax.set_xlabel(cvs[1])
    ax.set_ylabel("FES [kJ/mol]")

    # Save the figure
    plt.savefig("fes.png")  # Can use plt.savefig or fig.savefig
    plt.close()  # Good practice to close the figure after saving


if __name__ == "__main__":
    main()