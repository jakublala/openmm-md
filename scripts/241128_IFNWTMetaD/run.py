from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run(
        filepath: str, 
        system: str, 
        biasfactor: int, 
        sigma_cv1: str, 
        sigma_cv2: str, 
        grid_min_cv1: str, 
        grid_min_cv2: str, 
        grid_max_cv1: str, 
        grid_max_cv2: str,
        output_dir: str
        ):
        

    DATE = '241128'
    FILEPATH = filepath
    # OUTPUT_DIR = f'../../data/241010_FoldingUponBinding/output/{system}/{DATE}-MetaD'
    OUTPUT_DIR = output_dir
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 1000

    # 1. PLUMED CONFIG
    # CVs are cmap, d in that order
    config = {
        'type': 'metad',
        'temperature': TEMPERATURE,
        'stride': 500,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'metad.pace': 500,
        'cvs': ['cmap', 'd'],
        'metad.sigma': f'{sigma_cv1},{sigma_cv2}', # "0.04,0.01"
        'metad.height': 1.25, # 1/2 * kBT
        'metad.grid_min': f'{grid_min_cv1},{grid_min_cv2}',
        'metad.grid_max': f'{grid_max_cv1},{grid_max_cv2}',
        'metad.grid_bin': "200,200",
        'metad.biasfactor': biasfactor,
        'upper_wall.at': 5, # keep this at UW=5, we are primarily looking at BIASFACTOR now
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': None,
        'spot2_residues': None
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
        split_chains=False if ('A-synuclein' in FILEPATH) or ('CD28' in FILEPATH) else True, # HACK
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=2,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)