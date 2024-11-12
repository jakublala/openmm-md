import pandas as pd

def get_sigmas(directory, target, binder, run, cvs):
    assert len(cvs) == 2, "Only 2D FES are supported"
    return [get_sigma(directory, target, binder, run, cv) for cv in cvs]

def get_sigma(directory, target, binder, run, cv):
    # get the sigma from the kernels file
    kernels_file = f"{directory}/{target}_{binder}.kernels"
    df = pd.read_table(
        kernels_file,
        dtype=float,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["time", "cmap", "d", "sigma_cmap", "sigma_d", "height", "logweight"]
    )
    sigma = df[f"sigma_{cv}"].iloc[-1]
    return sigma.item()
