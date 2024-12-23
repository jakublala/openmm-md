from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run(
        filepath: str, 
        output_dir: str,
        metad_pace: int
        ):
        
    FILEPATH = filepath
    OUTPUT_DIR = output_dir
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    PADDING = 4
    MDTIME = 1000
    UPPER_WALL_AT = 5
    # 1. PLUMED CONFIG
    # CVs are cmap, d in that order
    config = {
        'type': 'metad',
        'temperature': TEMPERATURE,
        'stride': 500,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'metad.pace': metad_pace,
        'cvs': ['cmap', 'd'],
        'metad.sigma': "0.15,0.27", # "0.04,0.01"
        'metad.height': 1.25, # 1/2 * kBT
        'metad.grid_min': "0,0",
        'metad.grid_max': "45,7",
        'metad.grid_bin': "200,200",
        'metad.biasfactor': 48,
        'upper_wall.at': UPPER_WALL_AT, # keep this at UW=5, we are primarily looking at BIASFACTOR now
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': None,
        'spot2_residues': None,
        'restart': False,
        'trajectory_logging': True
    }

    import os
    from src.utils import get_gpu_indices
    if 'CUDA_VISIBLE_DEVICES' in os.environ:
        gpu_indices = get_gpu_indices()
    else:
        gpu_indices = None
    
    # 2. RUN MINIMIZATION AND SIMULATION
    main(
        filepath=FILEPATH,
        output_dir=OUTPUT_DIR,
        temperature=TEMPERATURE,
        mdtime=MDTIME,
        timestep=TIMESTEP,
        device_index=gpu_indices,
        device='cuda',
        split_chains=True,
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=PADDING,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)