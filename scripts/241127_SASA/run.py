from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run():
    DATE = '241202-SASA'
    FILEPATH = '../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb'
    OUTPUT_DIR = f'../../data/241010_FoldingUponBinding/output/CD28-G/{DATE}'
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 50
    PADDING = 4
    UPPER_WALL = 5

    # 1. PLUMED CONFIG
    config = {
        'type': 'opes-explore',
        'opes.pace': 500,
        'opes.barrier': 200,
        'temperature': TEMPERATURE,
        'stride': 500,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'upper_wall.at': UPPER_WALL,
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': None,
        'spot2_residues': None,
        'cvs': ['d', 'sasa'],
        'sasa.algo': 'HASEL',
        'sasa.spot_id': 2, # binder: 1, target: 2
        'sasa.stride': 10,
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
        padding=PADDING,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)