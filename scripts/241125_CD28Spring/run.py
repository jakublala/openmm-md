from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run():
    DATE = '241125-Spring'
    FILEPATH = '../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb'
    OUTPUT_DIR = f'../../data/241010_FoldingUponBinding/output/CD28-G/{DATE}'
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 500
    UPPER_WALL = 5
    PADDING = 4

    # 1. PLUMED CONFIG
    config = {
        'type': 'opes',
        'opes.pace': 500,
        'opes.barrier': 200,
        'temperature': TEMPERATURE,
        'stride': 500,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'upper_wall.at': UPPER_WALL,
        'upper_wall.exp': 6,
        'upper_wall.kappa': 250.0,
        'spot1_residues': None,
        'spot2_residues': None
    }
    
    # 1/4 of the upper wall at d=a+1
    spring_k = 2 * config['upper_wall.at'] / (config['upper_wall.at'] + 1)**2
    config['spring.k'] = spring_k * 0.25
    config['spring'] = True

    # 2. RUN MINIMIZATION AND SIMULATION
    main(
        filepath=FILEPATH,
        output_dir=OUTPUT_DIR,
        temperature=TEMPERATURE,
        mdtime=MDTIME,
        timestep=TIMESTEP,
        device_index='0,1',
        device='cuda',
        split_chains=False if ('A-synuclein' in FILEPATH) or ('CD28' in FILEPATH) else True, # HACK
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=PADDING,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)