import numpy as np
from typing import List
from tqdm import tqdm
import logging

logger = logging.getLogger(__name__)

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
        if not np.isclose(integral, 1.0, rtol=1e-2):
            logger.warning("Density is not normalized to 1, but it is %f instead", integral)
            density /= integral

        return density

    def _process_chunk(self, chunk_data):
        """Worker function for parallel KDE computation."""
        start, end = chunk_data
        data_chunk = self.data[start:end]
        weights_chunk = self.weights[start:end]
        
        diff_x = self._X[:, :, None] - data_chunk[:, 0]
        diff_y = self._Y[:, :, None] - data_chunk[:, 1]
        kernel = np.exp(-0.5 * ((diff_x / self.sigma[0])**2 + 
                            (diff_y / self.sigma[1])**2))
        return np.sum(weights_chunk * kernel, axis=2)
    

    def compute_kde_2d(
        self,
        grid_x: np.ndarray,
        grid_y: np.ndarray,
        chunk_size: int = 1000,
        n_processes: int = None  # Number of processes to use (None = use all available)
    ) -> np.ndarray:
        """Compute 2D kernel density estimation with parallel processing.
        
        Args:
            grid_x: Grid points for x dimension
            grid_y: Grid points for y dimension
            chunk_size: Number of data points to process at a time
            n_processes: Number of processes to use (defaults to CPU count)
        """
        from multiprocessing import Pool, cpu_count
        # Use fewer processes (half of CPU cores)
        n_processes = min(cpu_count() // 2, 16) if n_processes is None else n_processes
        
        density = np.zeros((len(grid_x), len(grid_y)))
        # Store meshgrid as instance variables for worker access
        self._X, self._Y = np.meshgrid(grid_x, grid_y, indexing='ij')
        
        # Precompute constant factors
        normalization = self.sigma[0] * self.sigma[1] * np.sum(self.weights) * 2 * np.pi
        
        # Prepare chunks for parallel processing
        num_points = self.data.shape[0]
        chunks = [(start, min(start + chunk_size, num_points)) 
                for start in range(0, num_points, chunk_size)]
        

        # Process chunks in parallel
        try:
            # Process chunks in parallel with timeout
            with Pool(n_processes) as pool:
                for partial_density in tqdm(
                    pool.imap(self._process_chunk, chunks),
                    total=len(chunks),
                    desc=f"Computing KDE ({n_processes} processes)"
                ):
                    density += partial_density
                    del partial_density
        except Exception as e:
            logger.error(f"KDE computation failed: {e}")
            raise
        finally:
            # Clean up
            del self._X, self._Y
        
        # Combine results
        density /= normalization
        
        # Check and fix normalization
        dx = grid_x[1] - grid_x[0]
        dy = grid_y[1] - grid_y[0]
        integral = np.sum(density) * dx * dy
        if not np.isclose(integral, 1.0, rtol=1e-2):
            logger.warning("Density is not normalized to 1, but it is %f instead", integral)
            density /= integral

        return density
