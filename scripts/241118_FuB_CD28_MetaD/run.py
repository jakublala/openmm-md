from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run(biasfactor: int, output_dir: str):
    FILEPATH = '../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb'
    OUTPUT_DIR = f'../../data/241010_FoldingUponBinding/output/CD28-G-MetaD/{output_dir}'
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 500

    # 1. PLUMED CONFIG
    # CVs are cmap, d in that order
    config = {
        'type': 'metad',
        'temperature': TEMPERATURE,
        'stride': 500,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'metad.pace': 100,
        'metad.sigma': "0.04,0.01",
        'metad.height': 0.1,
        'metad.grid_min': "0,0",
        'metad.grid_max': "80,5",
        'metad.grid_bin': "200,200",
        'metad.biasfactor': biasfactor,
        'upper_wall.at': 5,
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
        split_chains=True,
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=2,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)