import os
import logging
from typing import Optional, Literal
import shutil
from src.plumed.opes import run_plumed
from src.relax import minimize
from src.fixer import fixer
from src.plumed.io import create_plumed_input
from src.analysis.utils import get_file_by_extension
from src.plumed.utils import get_checkpoint_interval, get_pace_from_metad, get_last_checkpoint_timestep, process_hills_for_restart
            
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
        device_precision='double',
        output_dir=None,
        padding=None,
        split_chains=None,
        logging_frequency=100,
        box_size=None,
        config=None, # only contains opes stuff for now
        chain_mode: Optional[Literal['single-chain', 'two-chain']] = None,
        equilibrate_only=False,
        ):
    
    # TODO: the config should be saved somewhere as a JSON!!!
    # then very useful down the line for analysis
    # as parsing, for instance, plumed.dat is a mess



    print("SPLIT CHAINS ARE", split_chains)

    if chain_mode is None:
        raise ValueError('Chain mode is required')

    if output_dir is None:
        raise ValueError('Output directory is required')

    try:
        get_file_by_extension(output_dir, '.out')
    except FileNotFoundError:
        pass
    else:
        raise FileExistsError(f"Output directory {output_dir} and its .out file already exists, refusing to overwrite! Delete .out file if you believe it is safe to do so.")

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
    input_dir = os.path.dirname(filepath)

    if 'equilibrated' in filename:
        logger.info('Input PDB file is equilibrated, assuming it comes from a previous run...')
        logger.info('Copying fixed, equilibrated and solvated pdb files to output directory')
        fixed_pdb = get_file_by_extension(input_dir, '_fixed.pdb')
        equilibrated_pdb = get_file_by_extension(input_dir, '_equilibrated.pdb')
        solvated_pdb = get_file_by_extension(input_dir, '_solvated.pdb')
        # copy all to the output_dir
        for file in [fixed_pdb, equilibrated_pdb, solvated_pdb]:
            shutil.copy(file, output_dir)
        filename = filename.replace('_equilibrated', '')
        

    assert f'output/{filename}' not in os.listdir(), f"Folder output/{filename} already exists. It might overwrite existing data!"

    if 'opes' in config['type']:
        if config['restart']:
            raise NotImplementedError("Restarting from checkpoint in OPES not implemented yet")
            if config['restart_rfile'] is not None:
                assert os.path.exists(config['restart_rfile']), f"File {config['restart_rfile']} does not exist"
                logger.info("Found a restart_rfile for PLUMED, going to restart the OPES simulation from .state and .chk from OpenMM")
                restart = True
            else:
                restart = False
        else:
            restart_checkpoint = None
    else:
        if config['restart']:

            # move hills file to output dir
            # chop off lines that are after the checkpoint
            # get the checkpoint file
            # figure out the timestep
            # HILLS file numbers every deposited kernel
            # .out prints every 100 ps (i.e. 50,000 time steps)
            # we finished at (were WALLtimed) at 221,900 ps, i.e. 221.9 ns
            # how often do we get a checkpoint?
            # we save every 1 nanoseconds, according to get_checkpoint_interval
            plumed_file = get_file_by_extension(input_dir, 'plumed.dat')
            metad_pace = get_pace_from_metad(plumed_file)
            assert metad_pace == config['metad.pace'], 'Previous MetaD pace does not match the current one.'
            
            hills_file = get_file_by_extension(input_dir, '.hills')
            out_file = get_file_by_extension(input_dir, '.out')
            
            checkpoint_interval = get_checkpoint_interval(timestep)
            last_checkpoint_timestep = get_last_checkpoint_timestep(out_file, checkpoint_interval)
            num_hills_before_checkpoint = last_checkpoint_timestep // config['metad.pace']
            time_of_last_hill = int(num_hills_before_checkpoint * config['metad.pace'] * timestep * 0.001) # in ps 
            new_hills_file_lines = process_hills_for_restart(hills_file, time_of_last_hill)
            with open(f'{output_dir}/{filename}.hills', 'w') as f:
                for line in new_hills_file_lines:
                    f.write(line)
            logger.info("Restarting MetaD as requested...")
            restart_checkpoint = get_file_by_extension(input_dir, '.chk')

            # TODO: you should probably also copy the .plumed FILE!!!! not re-built it!!


        else:
            restart_checkpoint = None

    logger.info(f'==================== Running {filename} ====================')
    logger.info(f"Running with timestep {timestep} fs and mdtime {mdtime} ns")
    if config['type'] in ['opes', 'opes-explore']:
        logger.info(f"Energy barrier {config['opes.barrier']} kJ/mol for OPES")
        logger.info(f"Pace {config['opes.pace']} steps of depositing bias in OPES.")
    elif config['type'] == 'metad':
        logger.info(f"Bias factor {config['metad.biasfactor']} for MetaD")
    else:
        raise NotImplementedError(f"Type {config['type']} not implemented")
    logger.info(f"CV1: {config['cv1.type']} with {config['cv1.sigma']=}, {config['cv1.grid_min']=}, {config['cv1.grid_max']=}, {config['cv1.grid_bin']=}")
    logger.info(f"CV2: {config['cv2.type']} with {config['cv2.sigma']=}, {config['cv2.grid_min']=}, {config['cv2.grid_max']=}, {config['cv2.grid_bin']=}")

    if not os.path.exists(f"{output_dir}/{filename}_equilibrated.cif"):
        logger.info('No equilibrated cif file found, checking whether we need to run relaxation...')
        if not os.path.exists(f'{output_dir}/{filename}_solvated.cif'):

            if split_chains:
                if ('CD28' in filepath) or ('A-synuclein' in filepath):
                    raise Exception("CD28 and A-synuclein are already split, you might be doing something wrong...")

            logger.info('No solvated pdb file found, running solvation...')
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
                device_precision=device_precision,
                output_dir=output_dir,
                padding=padding,
                box_size=box_size,
                )
        else:
            logger.info('Solvated and equilibrated pdb files found, skipping solvation and relaxation')

        
    run_plumed(
        filename=filename, 
        mdtime=mdtime, 
        device_index=str(device_index),
        timestep=timestep,
        temperature=temperature,
        device=device,
        device_precision=device_precision,
        output_dir=output_dir,
        logging_frequency=logging_frequency,
        plumed_config=config,
        plumed_mode=chain_mode,
        restart_checkpoint=restart_checkpoint,
        equilibrate_only=equilibrate_only
        )



import fire
if __name__ == '__main__':
    fire.Fire(main)