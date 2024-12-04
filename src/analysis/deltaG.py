import numpy as np

def marginalize_fes(
        fes_2d: np.ndarray,
        kBT: float,
        dx: float,
        axis: int = 1
    ) -> np.ndarray:
    """
    Marginalize a 2D FES to 1D
    
    Parameters:
    fes_2d: 2D numpy array containing free energy values
    kBT: thermal energy (kB * T)
    dx: bin size of the CV
    axis: axis to marginalize over (i.e. the CV to keep)

    Returns:
    1D numpy array with marginalized free energy
    """
    # Convert FES to probability distribution
    prob_2d = np.exp(-fes_2d / kBT)
    
    # Integrate (sum) over the dimension to remove
    # the * dx accounts for the bin size, but is effectively unnecessary
    # as it leads to a constant factor shift in the FES, which gets removed lower
    prob_1d = np.sum(prob_2d * dx, axis=axis)
    
    # Convert back to free energy
    fes_1d = -kBT * np.log(prob_1d)
    
    # Shift minimum to zero if desired
    fes_1d = fes_1d - np.min(fes_1d)
    
    return fes_1d


def compute_deltaG(
        fes_1d: np.ndarray, 
        cv_bins: np.ndarray,
        kBT: float,
        threshold: float # in nm
    ) -> float:
    """
    Compute the free energy difference between the bound and unbound states determined by the threshold.
    
    Parameters:
    fes_1d: 1D free energy surface
    cv_bins: bin centers of the CV
    kBT: thermal energy (kB * T)
    threshold: distance threshold separating bound/unbound states (nm)
    
    Returns:
    deltaG: Free energy of binding (same units as kBT)
    """
    # Create masks for bound and unbound states
    bound_mask = cv_bins <= threshold
    unbound_mask = cv_bins > threshold
    
    # Convert FES to probabilities
    probabilities = np.exp(-fes_1d / kBT)
    
    # Calculate partition functions for bound and unbound states
    Z_bound = np.sum(probabilities[bound_mask])
    Z_unbound = np.sum(probabilities[unbound_mask])

    # Calculate deltaG = -kBT * ln(Z_bound/Z_unbound)
    deltaG = -kBT * np.log(Z_bound/Z_unbound)
    
    return deltaG