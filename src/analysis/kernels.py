import pandas as pd
from src.analysis.utils import get_file_by_extension

def get_sigmas(directory, target, binder, run, cvs):
    assert len(cvs) == 2, "Only 2D FES are supported"
    return [get_sigma(directory, target, binder, run, cv) for cv in cvs]

def get_sigma(directory, cv):
    try:
        kernels_file = get_file_by_extension(directory, '.kernels')
    except FileNotFoundError:
        kernels_file = get_file_by_extension(directory, '.hills')

    # get first line of kernels file
    with open(kernels_file, 'r') as file:
        labels = file.readline().split()[2:]

    df = pd.read_table(
        kernels_file,
        dtype=float,
        sep=r"\s+",
        comment="#",
        header=None,
        names=labels
    )
    sigma = df[f"sigma_{cv}"].iloc[-1]
    return sigma.item()
