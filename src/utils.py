from openmm.app import StateDataReporter
import os
import subprocess

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


def get_gpu_indices():
    """
    Convert GPU UUIDs from CUDA_VISIBLE_DEVICES to numerical indices for OpenMM.
    
    Returns
    -------
    str
        Comma-separated string of GPU indices (e.g., "0,1").
        Returns "0" if no GPUs are found or if there's an error.
    
    Examples
    --------
    >>> # If CUDA_VISIBLE_DEVICES="GPU-6aae2784-3658-0a4c-b8e4-510d2197db82,GPU-bebac720-9854-278d-fb9f-e3c0e55cad8e"
    >>> get_gpu_indices()
    '0,1'
    """
    def get_gpu_index_from_uuid(uuid):
        try:
            result = subprocess.run(['nvidia-smi', '-L'], capture_output=True, text=True)
            gpu_list = result.stdout.strip().split('\n')
            
            for i, gpu_info in enumerate(gpu_list):
                if uuid in gpu_info:
                    return str(i)
            return '0'
        except:
            return '0'

    try:
        # Get GPU UUIDs from environment
        gpu_uuids = os.environ.get('CUDA_VISIBLE_DEVICES', '').split(',')
        
        # Convert UUIDs to indices
        gpu_indices = []
        for uuid in gpu_uuids:
            uuid = uuid.strip()
            if uuid.startswith('GPU-'):
                index = get_gpu_index_from_uuid(uuid)
                gpu_indices.append(index)
        
        return ','.join(gpu_indices) if gpu_indices else '0'
    except:
        return '0'  # Default to GPU 0 if anything fails