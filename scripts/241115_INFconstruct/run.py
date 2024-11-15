import MDAnalysis as mda
from src.plumed.main import main
from src.plumed.io import create_plumed_input
from src.plumed.utils import get_checkpoint_interval
from src.models import Segment, Residue

if __name__ == '__main__':
    # Z-B50W
    FILEPATH = '../../data/241109_INFconstruct/input/Z1-B50W.pdb'
    OUTPUT_DIR = '../../data/241109_INFconstruct/output/Z1-B50W/241114'
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    MDTIME = 2

    # 1. CREATE PLUMED INPUT

    BINDER_LENGTH = 50
    LINKER1_LENGTH = 10
    PROTEASE_LENGTH = 10
    LINKER2_LENGTH = 10
    INF_LENGTH = 161

    NON_INF_LENGTH = BINDER_LENGTH + LINKER1_LENGTH + PROTEASE_LENGTH + LINKER2_LENGTH

    nresidues = mda.Universe(FILEPATH)._topology.n_residues
    assert nresidues == (
        NON_INF_LENGTH + INF_LENGTH
        ), f"Expected {NON_INF_LENGTH + INF_LENGTH} residues, got {nresidues}"


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
            NON_INF_LENGTH + 55 + 1,
            NON_INF_LENGTH + 135 + 1
            )
    ])

    config = {
        'type': 'opes-explore',
        'pace': 500,
        'barrier': 200,
        'temperature': TEMPERATURE,
        'stride': 500,
        'cutoff': 0.8,
        'restart_rfile': None,
        'state_wstride': get_checkpoint_interval(TIMESTEP),
        'upper_wall.at': 5,
        'upper_wall.exp': 6,
        'upper_wall.kappa': 1000.0,
        'spot1_residues': spot1_residues,
        'spot2_residues': spot2_residues
    }

    # 2. RUN MINIMIZATION AND SIMULATION
    main(
        filepath=FILEPATH,
        output_dir=OUTPUT_DIR,
        temperature=TEMPERATURE,
        mdtime=MDTIME,
        timestep=TIMESTEP,
        device_index='0,1',
        device='cuda',
        split_chains=False,
        logging_frequency=LOGGING_FREQUENCY,
        config=config,
        chain_mode='single-chain'
    )
