import numpy as np
import matplotlib.pyplot as plt

from src.analysis.utils import get_file_by_extension
from src.analysis.io import read_hills_file

def plot_hills(directory, ax):
    hills_file = get_file_by_extension(directory, '.hills')
    hills_df = read_hills_file(hills_file)
    ax.plot(hills_df['time'], hills_df['height'])
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Height [kJ/mol]")
    return ax