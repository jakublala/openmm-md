import numpy as np
import matplotlib.pyplot as plt

from src.analysis.utils import get_file_by_extension
from src.analysis.io import read_hills_file
from src.utils import get_latest_restart_dir

def plot_hills(directory, ax):
    directory = get_latest_restart_dir(directory)
    hills_file = get_file_by_extension(directory, '.hills')
    hills_df = read_hills_file(hills_file)
    ax.plot(hills_df['time'] / 1000, hills_df['height']) # convert ps to ns
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Height [kJ/mol]")
    return ax