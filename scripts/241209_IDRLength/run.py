from src.vanilla.main import main
from src.plumed.utils import get_checkpoint_interval
import fire

def run(
        filepath: str,
        output_dir: str,
        ):
        
    FILEPATH = filepath
    OUTPUT_DIR = output_dir
    TEMPERATURE = 300
    LOGGING_FREQUENCY = 100
    TIMESTEP = 2
    PADDING = 4
    MDTIME = 1000
    UPPER_WALL_AT = 5
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


    import os
    from src.utils import get_gpu_indices
    if 'CUDA_VISIBLE_DEVICES' in os.environ:
        gpu_indices = get_gpu_indices()
    else:
        gpu_indices = None
    
    # 2. RUN MINIMIZATION AND SIMULATION
    from src.vanilla.main import main
    main(
        filepath=FILEPATH,
        device_index=gpu_indices,
        output_dir=OUTPUT_DIR,
        temperature=TEMPERATURE,
        mdtime=MDTIME,
        timestep=TIMESTEP,
        device='cuda',
        split_chains=False,
        logging_frequency=LOGGING_FREQUENCY,
        padding=PADDING,
    )


if __name__ == '__main__':
    fire.Fire(run)