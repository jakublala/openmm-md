from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run():
    DATE = '241205-MetaDSASA'
    FILEPATH = '../../data/241010_FoldingUponBinding/output/CD28-G/241128-MetaD/CD28_general_equilibrated.pdb'
    OUTPUT_DIR = f'../../data/241010_FoldingUponBinding/output/CD28-G/{DATE}'
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 1000
    PADDING = 4
    UPPER_WALL = 5

    import MDAnalysis as mda
    universe = mda.Universe(FILEPATH)
    BINDER_LENGTH = len(universe.select_atoms('chainid A'))


    # FOR CD28
    # we will exclude this from computing the centre of masses
    if 'CD28' not in FILEPATH:
        raise NotImplementedError("Only CD28 is supported at the moment")
    from src.models import Segment, Residue
    idr_start_id = 119
    idr_end_id = 140
    idr_residues = Segment(residues=[
        Residue(
            index=i, 
            global_index=BINDER_LENGTH + i,
            chain_id='B',
            indexing=1
        ) for i in range(
            idr_start_id,
            idr_end_id + 1
            )
    ])

    # 1. PLUMED CONFIG
    config = {
        'type': 'metad',
        'metad.pace': 500,
        'metad.sigma': f'0.5,0.02',
        'metad.grid_min': '0,0',
        'metad.grid_max': f'60,7',
        'metad.grid_bin': '200,200',
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
        'idr_residues': idr_residues,
        'cvs': ['sasa', 'd'],
        'sasa.algo': 'HASEL',
        'sasa.spot_id': 2, # binder: 1, target: 2
        'sasa.stride': 10,
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
        split_chains=False if ('A-synuclein' in FILEPATH) or ('CD28' in FILEPATH) else True, # HACK
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=PADDING,
        chain_mode='two-chain'
    )


if __name__ == '__main__':
    fire.Fire(run)