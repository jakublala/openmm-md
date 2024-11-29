import pandas as pd
import re
import logging

logger = logging.getLogger(__name__)

def plot_colvar_trajectories(df, axs, timestep, stride):
    ax1, ax2 = axs

    # skip every 100th point, to smooth out the trajectories
    df = df.iloc[::100]

    # TODO: make this more general, cmap and d are hard-locked
    ax1.plot(df["time"] * timestep * stride, df["cmap"], label='CMAP')
    ax2.plot(df["time"] * timestep * stride, df["d"], label='Distance')

    ax1.set_title("CMAP")
    ax1.set_xlabel("Time [ns]")
    ax1.set_ylabel("CMAP [#]")
    ax1.legend()

    ax2.set_title("Distance")
    ax2.set_xlabel("Time [ns]")
    ax2.set_ylabel("Distance [nm]")
    ax2.legend()
    return ax1, ax2