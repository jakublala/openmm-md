import pandas as pd
import numpy as np
import h5py
import logging
from typing import List, Optional, Tuple
from src.constants import kB
from src.analysis.kde import GaussianKDE
# move this elsewhere then
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def compute_kde_weights(colvar_df: pd.DataFrame, bias: Optional[List[str]], kbT: float) -> np.ndarray:
    """Compute weights for KDE from bias columns"""
    total_bias = np.zeros(len(colvar_df))
    if bias is None:
        total_bias = colvar_df[[col for col in colvar_df.columns if 'bias' in col.lower()]].sum(axis=1)
    else:
        total_bias = colvar_df[bias].sum(axis=1)
    
    # Subtract maximum value to prevent overflow
    max_bias = np.max(total_bias.values)
    return np.exp((total_bias.values - max_bias) / kbT)

def compute_1d_fes(colvar_df: pd.DataFrame, cv: str, kde: GaussianKDE, n_bins: int, kbT: float) -> Tuple[np.ndarray, np.ndarray]:
    """Compute 1D FES using KDE"""
    cv_bins = np.linspace(colvar_df[cv].min(), colvar_df[cv].max(), n_bins)
    # Evaluate KDE for all points at once
    fes = -kbT * np.log(kde(cv_bins))
    fes -= np.min(fes)
    return cv_bins, fes

def compute_2d_fes(colvar_df: pd.DataFrame, cvs: List[str], kde: GaussianKDE, n_bins: int, kbT: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute 2D FES using KDE"""
    cv1_bins = np.linspace(colvar_df[cvs[0]].min(), colvar_df[cvs[0]].max(), n_bins)
    cv2_bins = np.linspace(colvar_df[cvs[1]].min(), colvar_df[cvs[1]].max(), n_bins)

    # Evaluate KDE for all points at once
    epsilon = 1e-300 # to avoid log(0)
    fes = -kbT * np.log(kde(cv1_bins, cv2_bins) + epsilon)
    fes = fes.reshape(n_bins, n_bins)
    fes -= np.min(fes)
    
    return cv1_bins, cv2_bins, fes

def save_fes(outfile: str, cv1_bins: np.ndarray, cv2_bins: Optional[np.ndarray], fes: np.ndarray, cvs: List[str]) -> None:
    """Save FES data to HDF5 file"""
    with h5py.File(outfile, 'w') as f:
        f.create_dataset(f'{cvs[0]}_bins', data=cv1_bins)
        if cv2_bins is not None:
            f.create_dataset(f'{cvs[1]}_bins', data=cv2_bins)
        f.create_dataset('fes', data=fes)
        f.attrs['description'] = 'Free Energy Surface'
        f.attrs['units'] = 'kJ/mol'

def compute_fes(
        colvar_df: pd.DataFrame,
        sigma: List[float],
        temp: float, 
        cvs: List[str],
        outfile: str,
        bias: Optional[List[str]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute free energy surface from collective variables and bias
    
    Args:
        colvar_df: DataFrame containing CVs and bias
        sigma: KDE bandwidth
        temp: Temperature in Kelvin
        cvs: List of CV column names
        outfile: Output HDF5 filename
        bias: Optional list of bias column names
        
    Returns:
        Tuple of (cv1_bins, cv2_bins, fes)
    """
    assert len(cvs) in [1, 2], "Only 1D and 2D FES are supported"
    assert len(sigma) == len(cvs), "Number of bandwidths must match number of CVs"

    N_BINS = 100
    kbT = kB * temp

    # Step 1: Compute the weights for the trajectory
    logger.info("Computing KDE weights")
    weights = compute_kde_weights(colvar_df, bias, kbT)

    # Step 2: Compute the KDE (unbiased probability density)
    logger.info("Computing the unbiased probability density")
    kde = GaussianKDE(colvar_df[cvs].values, weights=weights, sigma=sigma)
    
    # Step 3: Compute the FES
    if len(cvs) == 1:
        logger.info("Computing 1D FES")
        cv1_bins, fes = compute_1d_fes(colvar_df, cvs[0], kde, N_BINS, kbT)
        cv2_bins = None
    else:
        logger.info("Computing 2D FES")
        cv1_bins, cv2_bins, fes = compute_2d_fes(colvar_df, cvs, kde, N_BINS, kbT)

    # Save results
    logger.info(f"Saving FES as {outfile}")
    save_fes(outfile, cv1_bins, cv2_bins, fes, cvs)

    return cv1_bins, cv2_bins, fes