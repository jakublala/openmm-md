from src.plumed.main import main
from src.utils import get_checkpoint_interval
import fire
import re
import MDAnalysis as mda
from src.models import Residue, Segment

def extract_from_system(system):
    # need to use regex to extract for instance
    # Q7-B30L4 -> BINDER_LENGTH = 30, LINKER1_LENGTH = 4, LINKER2_LENGTH = 4
    # Z1-B40L10W -> BINDER_LENGTH = 40, LINKER1_LENGTH = 10, LINKER2_LENGTH = 10
    # PQ19-B50L10W -> BINDER_LENGTH = 50, LINKER1_LENGTH = 10, LINKER2_LENGTH = 10 
    pattern = r'([A-Z0-9]+)-B(\d+)L(\d+)W?'  # Made W optional with ? and added digits to first group
    match = re.match(pattern, system)
    if match:
        BINDER_LENGTH = int(match.group(2))
        LINKER1_LENGTH = int(match.group(3))
        LINKER2_LENGTH = LINKER1_LENGTH
    else:
        raise ValueError(f"System {system} does not match pattern {pattern}")
    return BINDER_LENGTH, LINKER1_LENGTH, LINKER2_LENGTH

def run(
        filepath: str, 
        system: str,
        output_dir: str,
        ):
        
    FILEPATH = filepath
    OUTPUT_DIR = output_dir
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 1000

    # 1. CREATE PLUMED INPUT
    if 'Q7' in system:
        PROTEASE_LENGTH = 9
    else:
        PROTEASE_LENGTH = 10
    BINDER_LENGTH, LINKER1_LENGTH, LINKER2_LENGTH = extract_from_system(system)
    INF_LENGTH = 161

    print(f"{LINKER1_LENGTH=}, {LINKER2_LENGTH=}, {BINDER_LENGTH=}, {system}")

    NON_INF_LENGTH = BINDER_LENGTH + LINKER1_LENGTH + PROTEASE_LENGTH + LINKER2_LENGTH

    # nresidues = mda.Universe(FILEPATH)._topology.n_residues
    from openmm.app import PDBxFile
    cif = PDBxFile(FILEPATH)
    universe = mda.Universe(cif)
    nresidues = len(universe.select_atoms("protein").residues)
    assert nresidues == (
        NON_INF_LENGTH + INF_LENGTH
        ), f"Expected {NON_INF_LENGTH + INF_LENGTH} residues, got {nresidues}"

    if 'W' in system:
        # this means a wider binding spot, 55 to 135
        spot2_start = 55
        spot2_end = 135
    else:
        # this means the default binding spot, 80 to 100
        spot2_start = 80
        spot2_end = 100

    # 1 - indexed
    spot1_residues = Segment(residues=[
        Residue(
            index=i, 
            global_index=i,
            chain_id='A', 
            indexing=1
        ) for i in range(1, BINDER_LENGTH + 1)
    ])
    spot2_residues = Segment(residues=[
        Residue(
            index=i, 
            global_index=i,
            chain_id='A',
            indexing=1
        ) for i in range(
            NON_INF_LENGTH + spot2_start + 1,
            NON_INF_LENGTH + spot2_end + 1
            )
    ])

    upper_wall_at = (3.8 * (LINKER1_LENGTH + PROTEASE_LENGTH + LINKER2_LENGTH)) / 10 / 2
    padding = upper_wall_at + 1

    # 1. PLUMED CONFIG
    # CVs are cmap, d in that order
    config = {
        'type': 'metad',
        'temperature': TEMPERATURE,
        'stride': 500,
        'cmap.contact_threshold': 0.8,
        'cmap.include_cutoff': 1.5,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'metad.pace': 500,
        'cv1.type': 'cmap',
        'cv1.sigma': 0.15,
        'cv1.grid_min': None,
        'cv1.grid_max': None,
        'cv1.grid_bin': None,
        'cv1.pbc': True,
        'cv2.type': 'd',
        'cv2.sigma': 0.27,
        'cv2.grid_min': None,
        'cv2.grid_max': None,
        'cv2.grid_bin': None,
        'cv2.pbc': True,
        'metad.height': 1.25, # 1/2 * kBT
        'metad.biasfactor': 50,
        'metad.height': 1.25, # 1/2 * kBT
        'upper_wall.at': upper_wall_at, # keep this at UW=5, we are primarily looking at BIASFACTOR now
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': spot1_residues,
        'spot2_residues': spot2_residues,
        'idr_residues': None,
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
        split_chains=False,
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        padding=padding,
        chain_mode='single-chain',
    )


if __name__ == '__main__':
    fire.Fire(run)