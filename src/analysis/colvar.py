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

    # remove the unbiased bit
    # PLUMED:   adaptive SIGMA will be used, with ADAPTIVE_SIGMA_STRIDE = 5000
    # PLUMED:     thus the first x kernel depositions will be skipped, x = ADAPTIVE_SIGMA_STRIDE/PACE = 10
    # get the pace from the _plumed.dat file
    filename = filename.replace(".colvar", "_plumed.dat")
    with open(filename, "r") as f:
        content = f.read()
    
    # depositing kernels
    pace_pattern = r"opes:\s*OPES_METAD[\s\S]*?PACE=(\d+)"
    match = re.search(pace_pattern, content)
    if match:
        pace = int(match.group(1))
        logger.info(f"Pace value: {pace}")
    else:
        raise ValueError("PACE not found in the file")
    
    # printing to COLVAR file
    stride_pattern = r"PRINT.*STRIDE=(\d+)"
    match = re.search(stride_pattern, content)
    if match:
        stride = int(match.group(1))
        logger.info(f"Stride value: {stride}")
    else:
        raise ValueError("Stride not found in the file")
    
    # assume first x = 10 kernel depositions are the ones to skip
    # compute how many printing steps correspond to that
    num_steps = pace * 10
    strides_to_skip = num_steps // stride
    logger.info(f"Skipping {strides_to_skip} steps")

    # look at the opes.bias in the strides_to_skip steps and assert they are constant
    assert colvar_df["opes.bias"].iloc[:strides_to_skip].std() < 1e-6, "opes.bias is not constant in non-biased part"
    colvar_df = colvar_df.iloc[strides_to_skip:]
    return colvar_df
