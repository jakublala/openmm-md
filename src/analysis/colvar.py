import pandas as pd
import re
import logging

logger = logging.getLogger(__name__)

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

    filename = filename.replace(".colvar", "_plumed.dat")
    with open(filename, "r") as f:
        content = f.read()
    
    # depositing kernels
    pace_pattern = r"PACE=(\d+)"
    matches = re.findall(pace_pattern, content)
    if not matches:
        raise ValueError("PACE not found in the file")
    if len(matches) > 1:
        raise ValueError(f"Multiple PACE values found in file: {matches}")
    
    pace = int(matches[0])
    logger.info(f"Pace value: {pace}")

    # printing to COLVAR file
    stride_pattern = r"PRINT.*STRIDE=(\d+)"
    match = re.search(stride_pattern, content)
    if match:
        stride = int(match.group(1))
        logger.info(f"Stride value: {stride}")
    else:
        raise ValueError("Stride not found in the file")
    
    # look at the opes.bias in the strides_to_skip steps and assert they are constant
    if 'opes.bias' in colvar_df.columns:
        # assume first x = 10 kernel depositions are the ones to skip
        # compute how many printing steps correspond to that
        num_steps = pace * 10
        strides_to_skip = num_steps // stride
        logger.info(f"Skipping {strides_to_skip} steps")
        assert colvar_df["opes.bias"].iloc[:strides_to_skip].std() < 1e-6, "opes.bias is not constant in non-biased part"
        colvar_df = colvar_df.iloc[strides_to_skip:]
    return colvar_df

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