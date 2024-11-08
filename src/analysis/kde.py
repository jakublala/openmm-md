import numpy as np
from typing import List
from tqdm import tqdm
class GaussianKDE:
    def __init__(self, data: np.ndarray, weights: np.ndarray, sigma: List[float]):
        self.data = data
        self.weights = weights
        self.sigma = np.array(sigma)
    
    def __call__(self, grid_x: np.ndarray, grid_y: np.ndarray = None) -> np.ndarray:
        if grid_y is None:
            return self.compute_kde_1d(grid_x)
        else:
            return self.compute_kde_2d(grid_x, grid_y)
        

    def compute_kde_1d(
        self,
        grid_x: np.ndarray
    ) -> np.ndarray:
        """Compute 1D kernel density estimation with a single bandwidth.
        
        Args:
            grid_x: Grid points for x dimension
        """
        """Compute 1D kernel density estimation with a single bandwidth.
        
        Args:
            grid_x: Grid points for x dimension
        
        Returns:
            density: 1D array of estimated density at grid points
        """
        # Compute distances and kernel for all points at once
        diff = grid_x[:] - self.data[:]
        kernel = np.exp(-0.5 * (diff / self.sigma[0])**2)
        density = np.sum(self.weights[:, None] * kernel, axis=0)
        
        # Normalize
        density /= self.sigma[0] * np.sum(self.weights) * np.sqrt(2 * np.pi)

        # Integrate over the grid points to check normalization
        dx = grid_x[1] - grid_x[0]
        integral = np.sum(density) * dx
        assert np.isclose(integral, 1.0, rtol=1e-2), f"Integral of density should be about 1, but it is {integral} instead"

        return density

    def compute_kde_2d(
        self,
        grid_x: np.ndarray,
        grid_y: np.ndarray,
        chunk_size: int = 10000
    ) -> np.ndarray:
        """Compute 2D kernel density estimation with different bandwidths for each dimension
        using chunking to reduce memory usage.
        
        Args:
            grid_x: Grid points for x dimension
            grid_y: Grid points for y dimension
            chunk_size: Number of data points to process at a time
            
        Returns:
            density: 2D array of estimated density at grid points
        """
        density = np.zeros((len(grid_x), len(grid_y)))
        
        # Create meshgrid of grid points
        X, Y = np.meshgrid(grid_x, grid_y, indexing='ij')  # Shape: (nx, ny)
        
        # Precompute constant factors
        normalization = self.sigma[0] * self.sigma[1] * np.sum(self.weights) * 2 * np.pi
        
        num_points = self.data.shape[0]
        for start in tqdm(range(0, num_points, chunk_size), desc="Computing KDE"):
            end = start + chunk_size
            data_chunk = self.data[start:end]
            weights_chunk = self.weights[start:end]
            
            # Compute differences
            diff_x = X[:, :, None] - data_chunk[:, 0]  # Shape: (nx, ny, chunk_size)
            diff_y = Y[:, :, None] - data_chunk[:, 1]  # Shape: (nx, ny, chunk_size)
            
            # Compute kernel
            kernel = np.exp(-0.5 * ((diff_x / self.sigma[0])**2 + (diff_y / self.sigma[1])**2))
            
            # Accumulate weighted kernel
            density += np.sum(weights_chunk * kernel, axis=2)
        
        # Normalize the density
        density /= normalization
        
        # Integrate over the grid points to check normalization
        dx = grid_x[1] - grid_x[0]
        dy = grid_y[1] - grid_y[0]
        integral = np.sum(density) * dx * dy
        assert np.isclose(integral, 1.0, rtol=1e-2), f"Integral of density should be about 1, but it is {integral} instead"
        
        return density
