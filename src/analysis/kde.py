import numpy as np
from typing import List

class GaussianKDE:
    def __init__(self, data: np.ndarray, weights: np.ndarray, sigma: List[float]):
        self.data = data
        self.weights = weights
        self.sigma = sigma
    
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
        grid_y: np.ndarray
    ) -> np.ndarray:
        """Compute 2D kernel density estimation with different bandwidths for each dimension.
        
        Args:
            data: Input data points (N x 2 array)
            grid_x: Grid points for x dimension
            grid_y: Grid points for y dimension
            sigma_x: Bandwidth for x dimension
            sigma_y: Bandwidth for y dimension
            weights: Weights for each data point
        
        Returns:
            density: 2D array of estimated density at grid points
        """
        # Create meshgrid for evaluation points
        X, Y = np.meshgrid(grid_x, grid_y)
        positions = np.vstack([X.ravel(), Y.ravel()])
        
        # Compute distances and kernel for all points at once
        diff = positions[:, None] - self.data[None, :]
        kernel = np.exp(-0.5 * (diff / self.sigma)**2)
        density = np.sum(self.weights[:, None] * kernel, axis=0)
        
        # Normalize
        density /= self.sigma[0] * self.sigma[1] * np.sum(self.weights) * np.sqrt(2 * np.pi)

        # Integrate over the grid points to check normalization
        dx = grid_x[1] - grid_x[0]
        dy = grid_y[1] - grid_y[0]
        integral = np.sum(density) * dx * dy
        assert np.isclose(integral, 1.0, rtol=1e-2), f"Integral of density should be about 1, but it is {integral} instead"

        return density
