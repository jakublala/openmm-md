import pandas as pd
import re
import logging
from io import StringIO

logger = logging.getLogger(__name__)

def read_colvar_file(filename, plumed_filename=None):
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

    if plumed_filename is None:
        plumed_filename = filename.replace(".colvar", "_plumed.dat")

    with open(plumed_filename, "r") as f:
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


def read_hills_file(filename):
    # Read first line to get column names
    with open(filename) as f:
        lines = f.readlines()
    
    # look if there's more than one line starting with #!
    header_lines = [line for line in lines if line.startswith('#!')]
    header_lines = [header_lines[i:i+3] for i in range(0, len(header_lines), 3)]

    if len(header_lines) > 1:
        # assert each element in header_lines is the same
        all_equal = all(header_lines[0] == header_line for header_line in header_lines)
        if not all_equal:
            raise ValueError("Header lines are not the same")
    
    # get the column names from all lines
    column_names = header_lines[0][0].replace('#! FIELDS', '').strip().split()
    logger.info(f"Column names: {column_names}")

    # Convert lines to string buffer, excluding header lines
    data_lines = [line for line in lines if not line.startswith('#!')]
    data_buffer = StringIO('\n'.join(data_lines))
    
    # Read data using column names from file
    hills_df = pd.read_table(
        data_buffer,
        dtype=float,
        sep=r"\s+",
        comment="#",
        header=None,
        names=column_names
    )

    if len(header_lines) > 1:
        # increment the timing of subsequent simulations
        # check whether the 'time' column is fully continous without any breaks
        # the difference between each row can be 1 or 100 or some other value
        if not len(hills_df['time'].diff().dropna().unique()) == 1:
            logger.info("Time column is not continuous, incrementing the simulation time")
            # find the difference between each row
            min_diff = hills_df['time'].diff().abs().dropna().unique().min()
            logger.info(f"Minimum difference (i.e. timestep of HILLS deposition): {min_diff}")
            
            # make all the time columns start at 0 and go up by min_diff
            hills_df['time'] = [i * min_diff for i in range(len(hills_df))]

            assert len(hills_df['time'].diff().dropna().unique()) == 1, "Time column is not continuous"

    # remove incomplete lines
    original_len = len(hills_df)
    hills_df = hills_df[hills_df.notna().all(axis=1)]
    removed_len = original_len - len(hills_df)
    if removed_len > 0:
        logger.info(f"Removed {removed_len} rows with NaN entries")
    return hills_df
