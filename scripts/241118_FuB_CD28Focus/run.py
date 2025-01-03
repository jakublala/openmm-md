from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run(system, padding, upper_wall):
    DATE = '241121'
    FILEPATH = '../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb'
    OUTPUT_DIR = f'../../data/241010_FoldingUponBinding/output/{system}/{DATE}'
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 500

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
        'upper_wall.at': upper_wall,
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': None,
        'spot2_residues': None
    }

    # 2. RUN MINIMIZATION AND SIMULATION
    main(
        filepath=FILEPATH,
        output_dir=OUTPUT_DIR,
        temperature=TEMPERATURE,
        mdtime=MDTIME,
        timestep=TIMESTEP,
        device_index='0',
        device='cuda',
        split_chains=False if ('A-synuclein' in FILEPATH) or ('CD28' in FILEPATH) else True, # HACK
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=padding,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)