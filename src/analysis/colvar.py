import pandas as pd
import re
import logging

logger = logging.getLogger(__name__)

def resolve_cv_label(cv):
    if cv == 'cmap':
        return 'CMAP [#]'
    elif cv == 'd':
        return 'Distance [nm]'
    elif cv == 'sasa':
        return 'SASA [nm^2]'
    else:
        raise ValueError(f"Unknown CV: {cv}")

def plot_colvar_trajectories(df, cvs, axs, timestep, stride):
    ax1, ax2 = axs

    # skip every 100th point, to smooth out the trajectories
    df = df.iloc[::100]

    # TODO: make this more general, cmap and d are hard-locked
    ax1.plot(df["time"] * timestep * stride, df[cvs[0]], label=cvs[0])
    ax2.plot(df["time"] * timestep * stride, df[cvs[1]], label=cvs[1])

    ax1.set_title(cvs[0])
    ax1.set_xlabel("Time [ns]")
    ax1.set_ylabel(resolve_cv_label(cvs[0]))
    ax1.legend()

    ax2.set_title(cvs[1])
    ax2.set_xlabel("Time [ns]")
    ax2.set_ylabel(resolve_cv_label(cvs[1]))
    ax2.legend()
    return ax1, ax2