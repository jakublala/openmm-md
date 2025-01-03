from openmm.app import StateDataReporter
import os
import subprocess
from openmm.unit import nanoseconds, picoseconds
from openmm import Platform
from openmm.app import PDBxFile, PDBFile
import logging
from src.analysis.utils import get_file_by_extension
logger = logging.getLogger(__name__)

def get_restarted_files_by_extension(directory, extension):
    restarted_files = [get_file_by_extension(directory, extension)] # first non-restarted original file
    restart_dirs = [f"{directory}/{d}" for d in os.listdir(directory) if d.startswith('restart-')]
    for d in restart_dirs:
        restarted_files.append(get_file_by_extension(d, extension))
    return restarted_files

def get_latest_restart_dir(directory):
    restart_dirs = [d for d in os.listdir(directory) if d.startswith('restart-')]
    if not restart_dirs:
        return directory
    return f"{directory}/restart-{max([int(d.split('-')[-1]) for d in restart_dirs])}"

def get_platform_and_properties(
        device, 
        device_index, 
        device_precision = 'mixed',
        ):
    """Get the platform and properties for the specified device."""
    if device == "cuda":
        logger.info(f'Using CUDA device {device_index}')
        platform = Platform.getPlatformByName('CUDA')
        if device_precision is None:
            properties = {'DeviceIndex': device_index}
        else:
            properties = {'DeviceIndex': device_index, 'Precision': device_precision}
    elif device == "cpu":
        logger.info('Using CPU')
        platform = Platform.getPlatformByName('CPU')
        properties = None
    elif device == "opencl":
        logger.info('Using OpenCL')
        platform = Platform.getPlatformByName('OpenCL')
        properties = {'Precision': 'mixed'}
    else:
        raise ValueError('Invalid device')
    logger.info(f'Platform used: {platform.getName()} with properties {properties}')
    return platform, properties

def get_checkpoint_interval(timestep):
    return int((1 * nanoseconds) / (timestep * 0.001 * picoseconds)) 
    # return int((0.001 * nanoseconds) / (timestep * 0.001 * picoseconds)) 

def get_full_reporter(filename, log_freq, nsteps):
    return StateDataReporter(
                file=f'tmp/{filename}.out',
                reportInterval=log_freq,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
                progress=True,
                remainingTime=True,
                totalSteps=nsteps
            )


def get_available_gpu_indices():
    """
    Get available GPU indices using nvidia-smi -L command.
    
    Returns
    -------
    str
        Comma-separated string of GPU indices (e.g., "0,1").
        Returns "0" if no GPUs are found or if there's an error.
    
    Examples
    --------
    >>> # If system has 2 GPUs
    >>> get_available_gpu_indices()
    '0,1'
    """
    try:
        result = subprocess.run(['nvidia-smi', '-L'], capture_output=True, text=True)
        if result.returncode != 0:
            return '0'
        
        # nvidia-smi -L output format is like:
        # GPU 0: NVIDIA GeForce RTX 2080 Ti (UUID: GPU-...)
        # GPU 1: NVIDIA GeForce RTX 2080 Ti (UUID: GPU-...)
        gpu_indices = []
        for line in result.stdout.strip().split('\n'):
            if line.startswith('GPU '):
                idx = line.split(':')[0].split(' ')[1]
                gpu_indices.append(idx)
        
        return ','.join(gpu_indices) if gpu_indices else '0'
    except (subprocess.SubprocessError, FileNotFoundError):
        return '0'

def save_equilibrated_state(
        simulation,
        output_dir,
        filename
        ) -> None:
    topology = simulation.topology
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBxFile.writeFile(topology, positions, open(f'{output_dir}/{filename}_equilibrated.cif', 'w'))
    logger.info(f'Equilibrated state saved to {output_dir}/{filename}_equilibrated.cif')


def convert_pdb_to_cif(output_dir: str, filename: str, file_type: str):
    """Convert PDB file to CIF format.
    
    Args:
        output_dir: Directory containing the files
        filename: Base name of the file without extension
        file_type: Type of file ('equilibrated' or 'solvated')
    """
    pdb_path = f"{output_dir}/{filename}_{file_type}.pdb"
    cif_path = f"{output_dir}/{filename}_{file_type}.cif"
    
    if os.path.exists(pdb_path):
        pdb = PDBFile(pdb_path)
        PDBxFile.writeFile(pdb.topology, pdb.positions, open(cif_path, 'w'))
        os.remove(pdb_path)  # Remove old PDB file


def get_gpu_utilization():
    """
    Get GPU utilization and memory usage for each GPU using nvidia-smi.
    
    Returns
    -------
    list[dict]
        List of dictionaries containing GPU info, one per GPU:
        [
            {
                'id': '0',
                'gpu_util': 45,  # GPU utilization percentage
                'memory_util': 80,  # Memory utilization percentage
                'memory_used': 8192,  # Memory used in MB
                'memory_total': 11019  # Total memory in MB
            },
            ...
        ]
    """
    try:
        # Format: index, gpu_util, memory_used, memory_total
        cmd = ['nvidia-smi', '--query-gpu=index,utilization.gpu,memory.used,memory.total', 
               '--format=csv,noheader,nounits']
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.warning("Failed to get GPU utilization")
            return []
        
        gpus = []
        for line in result.stdout.strip().split('\n'):
            # Parse CSV output
            idx, gpu_util, mem_used, mem_total = [x.strip() for x in line.split(',')]
            
            gpus.append({
                'id': idx,
                'gpu_util': int(gpu_util),
                'memory_used': int(mem_used),
                'memory_total': int(mem_total),
                'memory_util': round(int(mem_used) / int(mem_total) * 100, 1)
            })
            
        return gpus
        
    except (subprocess.SubprocessError, FileNotFoundError) as e:
        logger.warning(f"Error getting GPU utilization: {e}")
        return []