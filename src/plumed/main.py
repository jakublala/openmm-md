import os
import logging
from typing import Optional, Literal
import shutil
from mpi4py import MPI
from src.plumed.opes import run_plumed
from src.relax import minimize
from src.fixer import fixer
from src.plumed.io import create_plumed_input
from src.analysis.utils import get_file_by_extension
from src.plumed.utils import get_pace_from_metad, get_last_checkpoint_timestep, process_hills_for_restart
from src.utils import get_checkpoint_interval
from src.plumed.replica import run_replica_plumed
import sys

def _setup_mpi():
    IS_MPI_INITIALIZED = MPI.Is_initialized()
    if IS_MPI_INITIALIZED:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        n_procs = comm.Get_size()
        print("MPI initialized, running in multi process mode with n_procs = ", n_procs)
    else:
        rank = 0
        n_procs = 1
        print("No MPI initialized, running in single process mode")

    return IS_MPI_INITIALIZED, rank, n_procs

def _setup_logging(rank):
    # Configure logging
    if rank == 0:
        logger = logging.getLogger()
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            stream=sys.stdout
        )
        # Silence MDAnalysis
        logging.getLogger('MDAnalysis').setLevel(logging.WARNING)
        # Silence any other noisy loggers
        logging.getLogger('parmed').setLevel(logging.WARNING)
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        # Keep OpenMMTools debug logging
        logging.getLogger('openmmtools').setLevel(logging.DEBUG)
        logging.getLogger('openmmtools.multistate').setLevel(logging.DEBUG)
    else:
        class NoOpLogger:
            def __getattr__(self, name):
                def method(*args, **kwargs):
                    pass
                return method

        logger = NoOpLogger()
    return logger


def main(
        filepath=None, 
        device_index='0',
        mdtime=100, # in ns
        timestep=2,
        temperature=300,
        device='cuda',
        device_precision='mixed',
        output_dir=None,
        padding=None,
        split_chains=None,
        logging_frequency=100,
        box_size=None,
        config=None, # only contains opes stuff for now
        chain_mode: Optional[Literal['single-chain', 'two-chain']] = None,
        equilibrate_only=False,
        replica_exchange=False,
        replica_type: Optional[Literal['rest2']] = None,
        swap_time=None,
        temperatures=None,
        generate_plumed_input=True,
        ):


    IS_MPI_INITIALIZED, rank, n_procs = _setup_mpi()

    logger = _setup_logging(rank)

    # RUN this only on one process
    if chain_mode is None:
        raise ValueError('Chain mode is required')

    if output_dir is None:
        raise ValueError('Output directory is required')
    input_dir = os.path.dirname(filepath)

    # Deduce whether we are restarting from a previous run
    try:
        out_file = get_file_by_extension(input_dir, '.out')
        # Check if file has more than 2 lines
        with open(out_file, 'r') as f:
            line_count = sum(1 for line in f)
        if line_count > 2:
            logger.info("Found existing .out file with content, setting restart to True")
            config['restart'] = True
    except FileNotFoundError:
        config['restart'] = False

    # make a new subfolder for the restart
    # first look whether there is already a restart folder (which could be restart-i where i is the index of the restart)
    # get all folders with 'restart-' in the name
    if config['restart']:
        restart_folders = [f for f in os.listdir(input_dir) if f.startswith('restart-')]
        if restart_folders:
            logger.info(f"Found {len(restart_folders)} restart folders, going to add a new one...")
            # get the index of the last restart folder
            restart_folders.sort(key=lambda x: int(x.split('-')[-1]))
            last_restart_index = int(restart_folders[-1].split('-')[-1])
            # make a new restart folder with the next index
            new_restart_index = last_restart_index + 1

            filepath = os.path.join(input_dir, f'restart-{last_restart_index}', os.path.basename(filepath))
            # asser that the file exists
            assert os.path.exists(filepath), f"File {filepath} does not exist, cannot restart properly..."
            output_dir = os.path.join(input_dir, f'restart-{new_restart_index}')

        else:
            output_dir = os.path.join(input_dir, 'restart-1')


    try:
        get_file_by_extension(output_dir, '.out')
    except FileNotFoundError:
        pass
    else:
        # if the file is empty continue and overwrite it, otherwise raise an error
        if os.path.getsize(get_file_by_extension(output_dir, '.out')) == 0:
            logger.info("Output directory and its .out file already exists, but it is empty. Overwriting...")
        else:
            raise FileExistsError(f"Output directory {output_dir} and its .out file already exists, refusing to overwrite! Delete .out file if you believe it is safe to do so.")

    os.makedirs(output_dir, exist_ok=True)

    if config is None:
        raise ValueError('Config is required')
    
    if isinstance(device_index, tuple) or isinstance(device_index, list):
        # this might happen if we send in CUDA_VISIBIBLE_DEVICES, which get converted to a tuple/list
        device_index = ",".join(str(x) for x in device_index)
    
    if filepath is None:
        raise ValueError('Filepath is required')


    filename = os.path.basename(filepath).split('.')[0]
    input_dir = os.path.dirname(filepath)

    if 'equilibrated' in filename:
        filename = filename.replace('_equilibrated', '')
        if rank == 0:
            logger.info('Input PDB file is equilibrated, assuming it comes from a previous run...')
            logger.info('Copying fixed, equilibrated and solvated pdb files to output directory')
            fixed_pdb = get_file_by_extension(input_dir, '_fixed.pdb')
            equilibrated_pdb = get_file_by_extension(input_dir, '_equilibrated.cif')
            solvated_pdb = get_file_by_extension(input_dir, '_solvated.cif')
            # copy all to the output_dir
            for file in [fixed_pdb, equilibrated_pdb, solvated_pdb]:
                try:
                    shutil.copy(file, output_dir)
                except shutil.SameFileError:
                    logger.warning(f"File {file} is the same as the one in the output directory, skipping copy... You are likely doing something manually. Be careful!")
            

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
        if config['restart'] and rank == 0:

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

            # copy the PLUMED file, but modify it
            from src.plumed.utils import prepare_plumed_file_for_restart
            prepare_plumed_file_for_restart(plumed_file, output_dir, filename)
            
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

            if os.path.exists(f"{output_dir}/{filename}_equilibrated.pdb"):
                # no longer supported, need to convert into cif
                from openmm.app import PDBxFile, PDBFile
                pdb = PDBFile(f"{output_dir}/{filename}_equilibrated.pdb")
                PDBxFile.writeFile(pdb.topology, pdb.positions, open(f"{output_dir}/{filename}_equilibrated.cif", 'w'))
                os.remove(f"{output_dir}/{filename}_equilibrated.pdb")  # Optional: remove old PDB file
            # do the same for _solvated.pdb
            if os.path.exists(f"{output_dir}/{filename}_solvated.pdb"):
                pdb = PDBFile(f"{output_dir}/{filename}_solvated.pdb")
                PDBxFile.writeFile(pdb.topology, pdb.positions, open(f"{output_dir}/{filename}_solvated.cif", 'w'))
                os.remove(f"{output_dir}/{filename}_solvated.pdb")  # Optional: remove old PDB file


        else:
            restart_checkpoint = None

    logger.info(f'==================== Running {filename} ====================')
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Running with timestep {timestep} fs and mdtime {mdtime} ns")
    if config['type'] in ['opes', 'opes-explore']:
        logger.info(f"Energy barrier {config['opes.barrier']} kJ/mol for OPES")
        logger.info(f"Pace {config['opes.pace']} steps of depositing bias in OPES.")
    elif config['type'] == 'metad':
        logger.info(f"Bias factor {config['metad.biasfactor']} for MetaD")
    else:
        raise NotImplementedError(f"Type {config['type']} not implemented")

    logger.info(f"CV1: {config['cv1.type']} with {config['cv1.sigma']=}") #, {config['cv1.grid_min']=}, {config['cv1.grid_max']=}, {config['cv1.grid_bin']=}")
    logger.info(f"CV2: {config['cv2.type']} with {config['cv2.sigma']=}") #,  {config['cv2.grid_min']=}, {config['cv2.grid_max']=}, {config['cv2.grid_bin']=}")

    if not os.path.exists(f"{output_dir}/{filename}_equilibrated.cif"):
        logger.info('No equilibrated cif file found, checking whether we need to run relaxation...')
        if not os.path.exists(f'{output_dir}/{filename}_solvated.cif') and rank == 0:

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



    # wait for all processes to be ready
    MPI.COMM_WORLD.Barrier()

    if replica_exchange:
        assert IS_MPI_INITIALIZED, "Replica exchange requires MPI"
        if not os.path.exists(f"{output_dir}/{filename}_equilibrated.cif"):
            logger.info('No equilibrated cif file found, running equilibriation...')
            # HACK: this is very hacky!!!!
            if rank == 0:
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
                    equilibrate_only=True
                )

        # TODO: hacking this as I am trying to isolate the bug
        run_replica_plumed(
            filename=filename, 
            mdtime=mdtime, 
            timestep=timestep,
            swap_time=swap_time,
            temperatures=temperatures,
            device=device,
            device_precision=device_precision,
            output_dir=output_dir,
            logging_frequency=logging_frequency,
            plumed_config=config,
            chain_mode=chain_mode,
            replica_type=replica_type,
        )
    else:
        assert rank == 0, "Usual single replica run doesn't support MPI, we only support a single process run."
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
            equilibrate_only=equilibrate_only,
            generate_plumed_input=generate_plumed_input
            )




import fire
if __name__ == '__main__':
    fire.Fire(main)