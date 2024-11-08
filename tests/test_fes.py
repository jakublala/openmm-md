import pytest
import numpy as np
import pandas as pd
import h5py
import os
import sys
from numba import njit
sys.path.append('..')
from src.analysis.fes import compute_fes
from src.constants import kB

@njit
def simple_nvt(x0, dt, n_steps, gradient_fn, gamma=1.0, kT=1.0):
    n_steps = int(n_steps)
    traj_x = []
    # Initialize velocity randomly from Maxwell-Boltzmann distribution
    v = np.random.normal(0, np.sqrt(kT))
    
    x = x0
    for n in range(n_steps):
        # calculate the force
        f = -gradient_fn(x)
        
        # First half of velocity update
        v = v + (f * dt - gamma * v * dt) / 2
        
        # Add random force from Langevin thermostat
        noise = np.sqrt(2 * gamma * kT * dt) * np.random.normal(0, 1)
        v = v + noise / 2
        
        # Position update
        x = x + v * dt
        
        # Second half of updates
        f = -gradient_fn(x)
        v = v + (f * dt - gamma * v * dt) / 2
        v = v + noise / 2
        
        traj_x.append(x)
        
    return x, np.array(traj_x)

def test_compute_fes():
    @njit
    def energy_fn(x):
        return x**2
    @njit
    def gradient_fn(x):
        return 2*x
    
    @njit
    def bias_fn(x):
        return 1.0 * (x**2 - 2.0**2)**2
    @njit
    def bias_gradient_fn(x):
        return 4.0 * (x**2 - 2.0**2) * x
    
    @njit
    def biased_sim_gradient_fn(x):
        return gradient_fn(x) + bias_gradient_fn(x)
    
    kbT = 5.0
    SIGMA = 0.1

    _, traj_x = simple_nvt(x0=0.0, dt=0.001, n_steps=1e7, gradient_fn=biased_sim_gradient_fn, gamma=1.0, kT=kbT)

    # CASE 1: no bias consideration
    colvar_df = pd.DataFrame({'x': traj_x})
    x_bins, y_bins, fes = compute_fes(colvar_df, sigma=[SIGMA], temp=kbT/kB, cvs=['x'], outfile='fes.h5py', bias=None)

    energy_plus_bias = energy_fn(x_bins) + bias_fn(x_bins)
    energy_plus_bias -= np.min(energy_plus_bias)

    # take difference between FES and bias, only in range [-2,2]
    mask = (x_bins >= -2) & (x_bins <= 2)
    diff = np.abs(fes[mask] - energy_plus_bias[mask])
    
    assert np.average(diff) < 1.0, "FES should match energy_plus_bias up to a constant shift in range [-2,2]"
    assert len(x_bins) == 100, "x_bins should have 100 bins"
    assert y_bins is None, "y_bins should be None"
    assert len(fes) == 100, "FES should have 100 bins"

    # CASE 2: unbiased FES
    colvard_df = pd.DataFrame({'x': traj_x, 'custom_bias': bias_fn(traj_x)})
    x_bins, y_bins, fes = compute_fes(colvard_df, sigma=[SIGMA], temp=kbT/kB, cvs=['x'], outfile='fes.h5py', bias=['custom_bias'])

    # Only check difference in range [-2,2]
    diff = np.abs(fes[mask] - energy_fn(x_bins[mask]))
    
    assert np.average(diff) < 1.0, "FES should match energy_fn up to a constant shift in range [-2,2]"
    assert len(x_bins) == 100, "x_bins should have 100 bins"
    assert y_bins is None, "y_bins should be None"
    assert len(fes) == 100, "FES should have 100 bins"

    # open the file and check the data
    with h5py.File('fes.h5py', 'r') as f:
        assert 'x_bins' in f.keys()
        assert 'fes' in f.keys()
        fes_h5py = f['fes'][:]

    assert np.allclose(fes, fes_h5py), "FES should match"

    # clean up
    os.remove('fes.h5py')


if __name__ == "__main__":
    test_compute_fes()




