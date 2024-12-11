from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run(
        filepath: str = "../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb", 
        output_dir: str = 'test',
        ):
        
    FILEPATH = filepath
    OUTPUT_DIR = output_dir
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 1000
    
    # PADDING = 4
    BOX_SIZE = [14, 14, 14]

    # 1. PLUMED CONFIG
    # CVs are cmap, d in that order
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


    config = {
        'type': 'metad',
        'temperature': TEMPERATURE,
        'stride': 500,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'metad.pace': 500,
        'cv1.type': 'cmap',
        'cv1.sigma': 0.15,
        'cv1.grid_min': 0,
        'cv1.grid_max': 45,
        'cv1.grid_bin': 200,
        'cv1.pbc': False,
        'cv2.type': 'd',
        'cv2.sigma': 0.27,
        'cv2.grid_min': 0,
        'cv2.grid_max': 7,
        'cv2.grid_bin': 200,
        'cv2.pbc': False,
        'metad.height': 1.25, # 1/2 * kBT
        'metad.biasfactor': 48,
        'upper_wall.at': None, # keep this at UW=5, we are primarily looking at BIASFACTOR now
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': None,
        'spot2_residues': None,
        'idr_residues': idr_residues,
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
        device_precision='mixed',
        split_chains=False,
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        box_size=BOX_SIZE,
        chain_mode='two-chain',
        equilibrate_only=True,
    )


if __name__ == '__main__':
    fire.Fire(run)