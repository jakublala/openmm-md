import pandas as pd
import numpy as np
import h5py
import logging
from typing import List, Optional, Tuple
from src.constants import kB
from src.analysis.kde import GaussianKDE
from tqdm import tqdm

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
        f.attrs['cvs'] = cvs
        f.attrs['1st_axis'] = cvs[0]
        f.attrs['2nd_axis'] = cvs[1] if len(cvs) > 1 else None


def load_fes(filepath: str):
    with h5py.File(filepath, "r") as f:
        cvs = f.attrs['cvs']
        fes = f["fes"][:]
        cv1_bins = f[f"{cvs[0]}_bins"][:]
        cv2_bins = f[f"{cvs[1]}_bins"][:]
        first_axis = f.attrs['1st_axis']
        second_axis = f.attrs['2nd_axis']

    assert first_axis in cvs and second_axis in cvs, "First and second axis must be in CVs"
    assert first_axis != second_axis, "First and second axis must be different"
    assert first_axis == cvs[0] and second_axis == cvs[1], "First axis must be the first CV and second axis must be the second CV"
    assert fes.shape[0] == len(cv1_bins) and fes.shape[1] == len(cv2_bins), "FES must have the same shape as the CV bins"
    return cvs, fes, cv1_bins, cv2_bins


def compute_fes_from_hills(
        hills_df: pd.DataFrame,
        temp: float,
        cvs: List[str],
        outfile: str,
        biasfactor: float,
        n_bins: int = 200
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute FES from hills file.
    
    Uses the relationship: lim_(t->inf) F(s)= - (T + deltaT) / deltaT  *  V(s, t)
    
    Based on PLUMEd, BIASFACTOR = (T + deltaT) / T

    Doing this naively with a loop, taking about 1 minute to do a 100 ns run.
    """
    # Create grid for evaluation
    cv1_bins = np.linspace(hills_df[cvs[0]].min(), hills_df[cvs[0]].max(), n_bins)
    cv2_bins = np.linspace(hills_df[cvs[1]].min(), hills_df[cvs[1]].max(), n_bins)
    X, Y = np.meshgrid(cv1_bins, cv2_bins)
    
    # Initialize potential grid
    V = np.zeros((n_bins, n_bins))

    # Get parameters for each CV
    sigma1 = hills_df[f'sigma_{cvs[0]}'].iloc[0]  # assuming constant sigma
    sigma2 = hills_df[f'sigma_{cvs[1]}'].iloc[0]
    
    # Sum up all Gaussian contributions
    for _, hill in tqdm(hills_df.iterrows(), total=len(hills_df), desc='Summing up hills...'):
        # Calculate distances from hill center to all grid points
        d1 = (X - hill[cvs[0]]) ** 2 / (2 * sigma1 ** 2)
        d2 = (Y - hill[cvs[1]]) ** 2 / (2 * sigma2 ** 2)
        
        # Add Gaussian contribution
        V += hill['height'] * np.exp(-(d1 + d2))
    
    # Convert bias potential to free energy
    deltaT = temp * (biasfactor - 1)
    fes = -((temp + deltaT) / deltaT) * V
    
    # Shift minimum to zero
    fes -= np.min(fes)
    
    # Save results
    if outfile is not None:
        save_fes(outfile, cv1_bins, cv2_bins, fes, cvs)
    
    return cv1_bins, cv2_bins, fes

def compute_fes(
        colvar_df: pd.DataFrame,
        sigmas: List[float],
        temp: float, 
        cvs: List[str],
        outfile: str,
        bias: Optional[List[str]] = None,
        n_bins: int = 200
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute free energy surface from collective variables and bias.
    Reweights the trajectory. First, we compute the weights for the trajectory,
    based on the Boltzmann distribution. Then, we compute the kernel density estimation
    of the unbiased probability density. Last, we invert it to get the free energy surface.
    
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
    assert len(sigmas) == len(cvs), "Number of bandwidths must match number of CVs"

    kbT = kB * temp

    # Step 1: Compute the weights for the trajectory
    logger.info("Computing KDE weights")
    weights = compute_kde_weights(colvar_df, bias, kbT)

    # Step 2: Compute the KDE (unbiased probability density)
    logger.info("Computing the unbiased probability density")
    kde = GaussianKDE(colvar_df[cvs].values, weights=weights, sigmas=sigmas)
    
    # Step 3: Compute the FES
    if len(cvs) == 1:
        logger.info("Computing 1D FES")
        cv1_bins, fes = compute_1d_fes(colvar_df, cvs[0], kde, n_bins, kbT)
        cv2_bins = None
    else:
        logger.info("Computing 2D FES")
        cv1_bins, cv2_bins, fes = compute_2d_fes(colvar_df, cvs, kde, n_bins, kbT)

    # Save results
    if outfile is not None:
        logger.info(f"Saving FES as {outfile}")
        save_fes(outfile, cv1_bins, cv2_bins, fes, cvs)

    return cv1_bins, cv2_bins, fes