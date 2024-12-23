from src.plumed.main import main
from src.plumed.utils import get_checkpoint_interval
import fire
import numpy as np

def run(
        filepath: str = "results/CD28_general_equilibrated.cif", 
        output_dir: str = 'results',
        ):
        
    FILEPATH = filepath
    OUTPUT_DIR = output_dir
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 200

    T_MIN = 300
    T_MAX = 360
    N_REPLICAS = 4

    
    # PADDING = 2
    BOX_SIZE = [14, 14, 14]

    # 1. PLUMED CONFIG
    # CVs are cmap, d in that order
    import MDAnalysis as mda
    from openmm.app import PDBxFile
    pdb = PDBxFile(FILEPATH)
    universe = mda.Universe(pdb)
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
        'cmap.contact_threshold': 0.8,
        'cmap.include_cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'metad.pace': 500,
        'cv1.type': 'cmap',
        'cv1.sigma': 0.15,
        'cv1.grid_min': None,
        'cv1.grid_max': None,
        'cv1.grid_bin': None,
        # 'cv1.grid_min': 0,
        # 'cv1.grid_max': 45,
        # 'cv1.grid_bin': 200,
        'cv1.pbc': True,
        'cv2.type': 'd',
        'cv2.sigma': 0.27,
        'cv2.grid_min': None,
        'cv2.grid_max': None,
        'cv2.grid_bin': None,
        # 'cv2.grid_min': 0,
        # 'cv2.grid_max': 12,
        # 'cv2.grid_bin': 200,
        'cv2.pbc': True,
        'metad.height': 1.25, # 1/2 * kBT
        'metad.biasfactor': 48,
        'upper_wall.at': 6, # keep this at UW=5, we are primarily looking at BIASFACTOR now
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

    from openmm.unit import kelvin
    TEMPERATURES = [T_MIN + (T_MAX - T_MIN) * (np.exp(float(i) / float(N_REPLICAS-1)) - 1.0) / (np.e - 1.0) for i in range(N_REPLICAS)]
    TEMPERATURES = [T * kelvin for T in TEMPERATURES]
    
    # 2. RUN MINIMIZATION AND SIMULATION
    main(
        filepath=FILEPATH,
        output_dir=OUTPUT_DIR,
        temperature=TEMPERATURE,
        mdtime=MDTIME,
        timestep=TIMESTEP,
        device_index=None,
        device='cuda',
        device_precision='mixed',
        split_chains=False,
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        box_size=BOX_SIZE,
        chain_mode='two-chain',
        replica_exchange=True,
        swap_time=1, # in ps
        temperatures=TEMPERATURES,
    )


if __name__ == '__main__':
    fire.Fire(run)