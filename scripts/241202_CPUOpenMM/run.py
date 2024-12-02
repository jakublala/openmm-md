from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run():
    FILEPATH = '../../data/241010_FoldingUponBinding/output/CD28-G/241128-MetaD/CD28_general_equilibrated.pdb'
    OUTPUT_DIR = f'../../data/241010_FoldingUponBinding/output/CD28-G/241202-CPU'
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 10
    UPPER_WALL = 5
    PADDING = 2

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
        'metad.sigma': "0.11,0.27",
        'metad.height': 1.25,
        'metad.grid_min': "0,0",
        'metad.grid_max': "41,7",
        'metad.grid_bin': "200,200",
        'metad.biasfactor': 40,
        'upper_wall.at': UPPER_WALL,
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': None,
        'spot2_residues': None,
        'restart': False
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
        device='cpu',
        split_chains=True,
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=PADDING,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)