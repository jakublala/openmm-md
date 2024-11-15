import numpy as np
import pandas as pd
import h5py
import os
import sys
from numba import njit
sys.path.append('..')
from src.analysis.fes import compute_fes
from src.constants import kB
from typing import Callable, Tuple
import matplotlib.pyplot as plt


@njit
def simple_nvt(x0: np.ndarray, dt: float, n_steps: int, gradient_fn: Callable, gamma: float = 1.0, kT: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    n_steps = int(n_steps)
    ndim = x0.size
    # Pre-allocate trajectory array instead of using a list
    traj = np.zeros((n_steps, ndim))
    
    # Initialize velocities randomly from Maxwell-Boltzmann distribution
    v = np.random.normal(0, np.sqrt(kT), size=ndim)
    
    x = x0.copy()
    for n in range(n_steps):
        # calculate the forces
        f = -gradient_fn(x)
        
        # First half of velocity update
        v = v + (f * dt - gamma * v * dt) / 2
        
        # Add random forces from Langevin thermostat
        noise = np.sqrt(2 * gamma * kT * dt) * np.random.normal(0, 1, size=ndim)
        v = v + noise / 2
        
        # Position update
        x = x + v * dt
        
        # Second half of updates
        f = -gradient_fn(x)
        v = v + (f * dt - gamma * v * dt) / 2
        v = v + noise / 2
        
        traj[n, :] = x

    return x, traj

def test_compute_fes_1d():
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

    _, traj_x = simple_nvt(x0=np.array([0.0]), dt=0.001, n_steps=1e7, gradient_fn=biased_sim_gradient_fn, gamma=1.0, kT=kbT)

    # CASE 1: no bias consideration
    colvar_df = pd.DataFrame({'x': traj_x.squeeze()})
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
    colvard_df = pd.DataFrame({'x': traj_x.squeeze(), 'custom_bias': bias_fn(traj_x.squeeze())})
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

def plot_analytically(energy_fn, bias_fn, x_range, y_range):
        # Create comprehensive visualization of energy landscapes
    x_bins = np.linspace(x_range[0], x_range[1], 100)
    y_bins = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x_bins, y_bins, indexing='ij')

    # Calculate 2D surfaces
    energy_2d = np.array([[energy_fn(np.array([x, y])) for x in x_bins] for y in y_bins])
    bias_2d = np.array([[bias_fn(np.array([x, y])) for x in x_bins] for y in y_bins])
    total_2d = energy_2d + bias_2d

    # Create figure with subplots
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    # 1D slices along x (y=0) with forces
    x_slice = np.array([energy_fn(np.array([x, 0])) for x in x_bins])
    bias_slice = np.array([bias_fn(np.array([x, 0])) for x in x_bins])
    total_slice = x_slice + bias_slice
    
    # Calculate forces (negative gradient)
    force_x = -np.gradient(x_slice, x_bins)
    force_bias = -np.gradient(bias_slice, x_bins)
    force_total = -np.gradient(total_slice, x_bins)

    # Plot energies
    ax1 = axs[0, 0]
    ax1.plot(x_bins, x_slice, 'b-', label='Energy')
    ax1.plot(x_bins, bias_slice, 'r-', label='Bias')
    ax1.plot(x_bins, total_slice, 'k--', label='Total')
    ax1.set_title('1D Slice (y=0)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('Energy')
    ax1.legend()

    # Add force subplot
    ax1_force = ax1.twinx()
    ax1_force.plot(x_bins, force_x, 'b:', alpha=0.5, label='Force (Energy)')
    ax1_force.plot(x_bins, force_bias, 'r:', alpha=0.5, label='Force (Bias)')
    ax1_force.plot(x_bins, force_total, 'k:', alpha=0.5, label='Force (Total)')
    ax1_force.set_ylabel('Force')
    ax1_force.legend(loc='center right')

    # 2D contour plots
    im0 = axs[0, 1].contourf(X, Y, energy_2d, levels=20)
    axs[0, 1].set_title('Energy Surface')
    axs[0, 1].set_xlabel('x')
    axs[0, 1].set_ylabel('y')
    axs[0, 1].set_aspect('equal')
    plt.colorbar(im0, ax=axs[0, 1])
    
    im1 = axs[0, 2].contourf(X, Y, bias_2d, levels=20)
    axs[0, 2].set_title('Bias Surface')
    axs[0, 2].set_xlabel('x')
    axs[0, 2].set_ylabel('y')
    axs[0, 2].set_aspect('equal')
    plt.colorbar(im1, ax=axs[0, 2])
    
    im2 = axs[1, 1].contourf(X, Y, total_2d, levels=20)
    axs[1, 1].set_title('Total Surface')
    axs[1, 1].set_xlabel('x')
    axs[1, 1].set_ylabel('y')
    axs[1, 1].set_aspect('equal')
    plt.colorbar(im2, ax=axs[1, 1])

    # Diagonal slice with forces
    diag_coords = np.array([[x, x] for x in x_bins])
    diag_energy = np.array([energy_fn(coord) for coord in diag_coords])
    diag_bias = np.array([bias_fn(coord) for coord in diag_coords])
    diag_total = diag_energy + diag_bias

    # Calculate forces along diagonal
    force_diag_energy = -np.gradient(diag_energy, x_bins)
    force_diag_bias = -np.gradient(diag_bias, x_bins)
    force_diag_total = -np.gradient(diag_total, x_bins)

    ax2 = axs[1, 0]
    ax2.plot(x_bins, diag_energy, 'b-', label='Energy')
    ax2.plot(x_bins, diag_bias, 'r-', label='Bias')
    ax2.plot(x_bins, diag_total, 'k--', label='Total')
    ax2.set_title('1D Slice (y=x)')
    ax2.set_xlabel('x=y')
    ax2.set_ylabel('Energy')
    ax2.legend()

    # Add force subplot for diagonal
    ax2_force = ax2.twinx()
    ax2_force.plot(x_bins, force_diag_energy, 'b:', alpha=0.5, label='Force (Energy)')
    ax2_force.plot(x_bins, force_diag_bias, 'r:', alpha=0.5, label='Force (Bias)')
    ax2_force.plot(x_bins, force_diag_total, 'k:', alpha=0.5, label='Force (Total)')
    ax2_force.set_ylabel('Force')
    ax2_force.legend(loc='center right')

    # Add circular mask visualization
    circle = plt.Circle((0, 0), 4, fill=False, color='red', linestyle='--')
    axs[1, 2].contourf(X, Y, total_2d, levels=20)
    axs[1, 2].add_patch(circle)
    axs[1, 2].set_title('Analysis Region (r≤4)')
    axs[1, 2].set_xlabel('x')
    axs[1, 2].set_ylabel('y')
    axs[1, 2].set_aspect('equal')
    plt.colorbar(im2, ax=axs[1, 2])

    plt.savefig('energy_landscapes.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_analyticall_with_fes(X, Y, fes, energy, energy_plus_bias, x_bins, y_bins):
    # Create figure with subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    # Calculate levels based on the analytical energy
    max_energy = np.max(energy)
    levels = np.linspace(0, max_energy, 20)

    # 2D contour plots
    im0 = axs[0, 0].contourf(X, Y, fes, levels=levels)
    axs[0, 0].set_title('Computed FES')
    plt.colorbar(im0, ax=axs[0, 0])
    axs[0, 0].set_aspect('equal')
    axs[0, 0].set_xlim(-3, 3)
    axs[0, 0].set_ylim(-3, 3)
    
    im1 = axs[0, 1].contourf(X, Y, energy, levels=levels)
    axs[0, 1].set_title('Analytical Energy')
    plt.colorbar(im1, ax=axs[0, 1])
    axs[0, 1].set_aspect('equal')
    axs[0, 1].set_xlim(-3, 3)
    axs[0, 1].set_ylim(-3, 3)

    # 1D slices
    # Slice along x (y=0)
    middle_idx = len(y_bins) // 2
    axs[1, 0].plot(x_bins, fes[middle_idx, :], 'b-', label='FES')
    axs[1, 0].plot(x_bins, energy[middle_idx, :], 'r--', label='Analytical')
    axs[1, 0].plot(x_bins, energy_plus_bias[middle_idx, :], 'g-.', label='Energy + Bias')
    axs[1, 0].set_title('1D Slice (y=0)')
    axs[1, 0].set_xlabel('x')
    axs[1, 0].set_ylabel('Energy')
    axs[1, 0].legend()
    axs[1, 0].set_xlim(-3, 3)
    axs[1, 0].set_ylim(-10, 40)

    # Diagonal slice
    diagonal = np.diagonal(fes)
    analytical_diagonal = np.diagonal(energy)
    analytical_diagonal_plus_bias = np.diagonal(energy_plus_bias)
    diag_coords = np.linspace(x_bins[0], x_bins[-1], len(diagonal))
    axs[1, 1].plot(diag_coords, diagonal, 'b-', label='FES')
    axs[1, 1].plot(diag_coords, analytical_diagonal, 'r--', label='Analytical')
    axs[1, 1].plot(diag_coords, analytical_diagonal_plus_bias, 'g-.', label='Energy + Bias')
    axs[1, 1].set_title('Diagonal Slice (y=x)')
    axs[1, 1].set_xlabel('x=y')
    axs[1, 1].set_ylabel('Energy')
    axs[1, 1].legend()
    axs[1, 1].set_xlim(-3, 3)
    axs[1, 1].set_ylim(-10, 40)

    plt.savefig('fes_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_trajectory(traj_x):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    
    # Time array (assuming uniform timesteps)
    time = np.arange(len(traj_x))
    
    # Plot x coordinate vs time
    ax1.scatter(time, traj_x[:, 0], s=0.1)
    ax1.set_ylabel('x')
    
    # Plot y coordinate vs time
    ax2.scatter(time, traj_x[:, 1], s=0.1)
    ax2.set_ylabel('y')
    ax2.set_xlabel('Time step')
    
    plt.tight_layout()
    plt.savefig('trajectory.png', dpi=300, bbox_inches='tight')
    plt.close()

def test_compute_fes_2d():
    @njit
    def energy_fn(x):
        return 0.25 * x[..., 0]**2 + 0.25 * x[..., 1]**2
    
    @njit
    def gradient_fn(x):
        grad = np.zeros_like(x)
        grad[0] = 0.5 * x[0]
        grad[1] = 0.5 * x[1]
        return grad
    
    @njit
    def bias_fn(x):
        # Double-well potential with minima at (2,2) and (-2,-2)
        # Using negative exponentials to create wells
        well1 = -5.0 * np.exp(-(((x[..., 0] - 2)**2 + (x[..., 1] - 2)**2)/1.0))
        well2 = -5.0 * np.exp(-(((x[..., 0] + 2)**2 + (x[..., 1] + 2)**2)/1.0))
        return well1 + well2
    
    @njit
    def bias_gradient_fn(x):
        # Gradient of the double-well potential
        grad = np.zeros_like(x)
        
        # Derivatives of exponential wells
        exp1 = np.exp(-(((x[0] - 2)**2 + (x[1] - 2)**2)/1.0))
        exp2 = np.exp(-(((x[0] + 2)**2 + (x[1] + 2)**2)/1.0))
        
        # dx component
        grad[0] = 10.0 * (x[0] - 2) * exp1 + 10.0 * (x[0] + 2) * exp2
        
        # dy component
        grad[1] = 10.0 * (x[1] - 2) * exp1 + 10.0 * (x[1] - 2) * exp2
        
        return grad
    
    @njit
    def biased_sim_gradient_fn(x):
        return gradient_fn(x) + bias_gradient_fn(x)
    
    if DEBUG:
        plot_analytically(energy_fn, bias_fn, [-6, 6], [-6, 6])

    kbT = 1.0
    SIGMA = [0.1, 0.1]

    _, traj_x = simple_nvt(x0=np.array([0.0, 0.0]), dt=0.001, n_steps=1e8, 
                                       gradient_fn=biased_sim_gradient_fn, gamma=1.0, kT=kbT)
    
    # only take 10th step of the trajectory
    traj_x = traj_x[::100, :]

    if DEBUG:
        plot_trajectory(traj_x)

    # CASE 1: no bias consideration
    colvar_df = pd.DataFrame({'x': traj_x[:, 0], 'y': traj_x[:, 1]})
    x_bins, y_bins, fes = compute_fes(colvar_df, sigma=SIGMA, temp=kbT/kB, 
                                     cvs=['x', 'y'], outfile='fes.h5py', bias=None)
    
    # Create meshgrid for evaluation
    X, Y = np.meshgrid(x_bins, y_bins, indexing='ij')
    coords = np.vstack([X.ravel(), Y.ravel()]).T
    energy_plus_bias = (energy_fn(coords) + bias_fn(coords)).reshape(len(x_bins), len(y_bins))
    energy_plus_bias -= np.min(energy_plus_bias)

    # take difference between FES and bias, only in circle of radius 3
    R = np.sqrt(X**2 + Y**2)
    mask = R <= 3
    diff = np.abs(fes - energy_plus_bias)[mask]
    
    assert np.average(diff) < 5.0, "FES should match energy_plus_bias up to a constant shift within radius 3"
    assert len(x_bins) == 100 and len(y_bins) == 100, "x_bins and y_bins should have 100 bins each"
    assert fes.shape == (100, 100), "FES should be 100x100"

    # CASE 2: unbiased FES
    bias_values = np.array([bias_fn(x) for x in traj_x])
    colvar_df = pd.DataFrame({
        'x': traj_x[:, 0], 
        'y': traj_x[:, 1],
        'custom_bias': bias_values
    })
    x_bins, y_bins, fes = compute_fes(colvar_df, sigma=SIGMA, temp=kbT/kB, 
                                     cvs=['x', 'y'], outfile='fes.h5py', 
                                     bias=['custom_bias'])

    energy = energy_fn(coords).reshape(len(x_bins), len(y_bins))

    if DEBUG:
        plot_analyticall_with_fes(X, Y, fes, energy, energy_plus_bias, x_bins, y_bins)

    # Only check difference within radius 2
    diff = np.abs(fes - energy)[mask]

    assert np.average(diff) < 5.0, "FES should match energy_fn up to a constant shift within radius 2"
    assert len(x_bins) == 100 and len(y_bins) == 100, "x_bins and y_bins should have 100 bins each"
    assert fes.shape == (100, 100), "FES should be 100x100"

    # Check the saved file
    with h5py.File('fes.h5py', 'r') as f:
        assert 'x_bins' in f.keys()
        assert 'y_bins' in f.keys()
        assert 'fes' in f.keys()
        fes_h5py = f['fes'][:]

    assert np.allclose(fes, fes_h5py), "FES should match"

    # clean up
    os.remove('fes.h5py')


if __name__ == "__main__":
    DEBUG = True
    # test_compute_fes_1d()
    test_compute_fes_2d()



