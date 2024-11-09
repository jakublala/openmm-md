import numpy as np

def marginalize_fes(
        fes_2d: np.ndarray,
        kBT: float,
        axis: int = 1
    ) -> np.ndarray:
    """
    Marginalize a 2D FES to 1D
    
    Parameters:
    fes_2d: 2D numpy array containing free energy values
    kBT: thermal energy (kB * T)
    axis: axis to marginalize over (i.e. the CV to keep)

    Returns:
    1D numpy array with marginalized free energy
    """
    # Convert FES to probability distribution
    prob_2d = np.exp(-fes_2d / kBT)
    
    # Integrate (sum) over the dimension to remove
    prob_1d = np.sum(prob_2d, axis=axis)
    
    # Convert back to free energy
    fes_1d = -kBT * np.log(prob_1d)
    
    # Shift minimum to zero if desired
    fes_1d = fes_1d - np.min(fes_1d)
    
    return fes_1d