import os
import logging
from typing import Optional, Literal

from src.plumed.opes import opes
from src.relax import minimize
from src.fixer import fixer
from src.plumed.io import create_plumed_input

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
logger = logging.getLogger(__name__)

def main(
        filepath=None, 
        device_index='0',
        mdtime=100, # in ns
        timestep=2,
        temperature=300,
        device='cuda',
        output_dir=None,
        padding=None,
        split_chains=None,
        logging_frequency=100,
        config=None, # only contains opes stuff for now
        chain_mode: Optional[Literal['single-chain', 'two-chain']] = None,
        ):
    
    if chain_mode is None:
        raise ValueError('Chain mode is required')

    if output_dir is None:
        raise ValueError('Output directory is required')
    os.makedirs(output_dir, exist_ok=True)

    if config is None:
        raise ValueError('Config is required')
    
    logger.info(f"Output directory: {output_dir}")

    if isinstance(device_index, tuple) or isinstance(device_index, list):
        # this might happen if we send in CUDA_VISIBIBLE_DEVICES, which get converted to a tuple/list
        device_index = ",".join(str(x) for x in device_index)
    
    if filepath is None:
        raise ValueError('Filepath is required')

    filename = os.path.basename(filepath).split('.')[0]

    assert f'output/{filename}' not in os.listdir(), f"Folder output/{filename} already exists. It might overwrite existing data!"

    if config['restart_rfile'] is not None:
        assert os.path.exists(config['restart_rfile']), f"File {config['restart_rfile']} does not exist"
        logger.info("Found a restart_rfile for PLUMED, going to restart the OPES simulation from .state and .chk from OpenMM")
        restart = True
    else:
        restart = False

    logger.info(f'==================== Running {filename} ====================')
    logger.info(f"Running with timestep {timestep} fs and mdtime {mdtime} ns")
    logger.info(f"Energy barrier {config['barrier']} kJ/mol for OPES")
    logger.info(f"Pace {config['pace']} steps of depositing bias in OPES.")

    
    if not restart:
        restart_checkpoint = None
        if not os.path.exists(f'{output_dir}/{filename}_solvated.pdb'):
            logger.info('No minimized pdb file found, running relaxation...')
            # 1. load the PDB and fix errors
            fixer(
                filepath=filepath, 
                output_dir=output_dir,
                split_chains=split_chains
                )
            logger.info("Fixing successful.")
            
            # 2. minimize the structure with LBFGS and H atoms mobile
            minimize(
                filename=filename, 
                max_iterations=0, 
                device_index=str(device_index),
                constraints=None,
                device=device,
                output_dir=output_dir,
                padding=padding
                )
            
    else:
        raise NotImplementedError("Restarting from checkpoint not implemented yet")
    
    create_plumed_input(
        filepath=filepath, 
        output_dir=output_dir,
        config=config,
        mode=chain_mode
        )


        
    opes(
        filename=filename, 
        mdtime=mdtime, 
        device_index=str(device_index),
        timestep=timestep,
        temperature=temperature,
        restart_checkpoint=restart_checkpoint,
        device=device,
        output_dir=output_dir,
        logging_frequency=logging_frequency
        )



import fire
if __name__ == '__main__':
    fire.Fire(main)