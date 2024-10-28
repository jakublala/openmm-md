import sys
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
import os
import subprocess
import shutil

logger = logging.getLogger(__name__)

# current datetime
from datetime import datetime

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error running command: {command}")
        print(f"Return code: {process.returncode}")
        print(f"Standard output:\n{stdout.decode()}")
        print(f"Standard error:\n{stderr.decode()}")
        raise Exception(f"Error running command: {command}")

    
    return stdout.decode()

def main(
        filepath=None, 
        device_index=1,
        mdtime=100, # in ns
        timestep=2,
        barrier=100,
        restart_rfile=None,
        pace=500,
        device='cuda',
        output_dir=None,
        ):

    if isinstance(device_index, tuple) or isinstance(device_index, list):
        # this might happen if we send in CUDA_VISIBIBLE_DEVICES, which get converted to a tuple/list
        device_index = ",".join(str(x) for x in device_index)
    
    CUTOFF = 1.0

    if filepath is None:
        raise ValueError('Filepath is required')

    filename = os.path.basename(filepath).split('.')[0]

    assert f'output/{filename}' not in os.listdir(), f"Folder output/{filename} already exists. It might overwrite existing data!"

    if restart_rfile is not None:
        assert os.path.exists(restart_rfile), f"File {restart_rfile} does not exist"
        logger.info(f"Found a restart_rfile for PLUMED, going to restart the OPES simulation from .state and previous last frame in trajectory.")
        restart = True
    else:
        restart = False

    logger.info(f'==================== Running {filename} ====================')
    if output_dir is None:
        logger.warning(f"WARNING: No output directory specified, using default tmp/{filename}")
        # check that it's empty
        if os.listdir(f'tmp/{filename}'):
            raise ValueError(f"Folder tmp/{filename} is not empty, please specify an output directory")
        else:
            os.makedirs(f'tmp/{filename}', exist_ok=True)
    else:
        os.makedirs(output_dir, exist_ok=True)

    if not restart:
        restart_checkpoint = None
        if not os.path.exists(f'{output_dir}/{filename}_solvated.pdb'):
            logger.info('No minimized pdb file found, running relaxation...')
            # 1. load the PDB and fix errors
            from src.fixer import fixer
            fixer(filepath=filepath, output_dir=output_dir)
            logger.info("Fixing successful.")
            

            # 2. minimize the structure with LBFGS and H atoms mobile
            from src.relax import minimize
            minimize(
                filename=filename, 
                max_iterations=0, 
                device_index=str(device_index),
                constraints=None,
                device=device,
                output_dir=output_dir
                )
            
        # 2.5 get OPES preparation
        from src.plumed.cv import get_interface_contact_indices
        contact_indices = get_interface_contact_indices(filename=filename, cutoff=CUTOFF, output_dir=output_dir)

        contact_pairs_str = ""
        for i, pair in enumerate(contact_indices):
            if i != 0:
                contact_pairs_str += f"\n\tATOMS{i+1}={pair[0]},{pair[1]}"
            else:
                contact_pairs_str += f"\tATOMS{i+1}={pair[0]},{pair[1]}"
    else:
        print(f"Getting CVs from previous plumed.dat")
        # assume it's in the same folder as restart_rfile
        restart_rfile_path = os.path.dirname(restart_rfile)
        plumed_file = os.path.join(restart_rfile_path, f'{filename}_plumed.dat')
        assert os.path.exists(plumed_file), f"File {plumed_file} does not exist"

        # assume fixed.pdb also in the same folder
        fixed_pdb = os.path.join(restart_rfile_path, f'{filename}_fixed.pdb')
        assert os.path.exists(fixed_pdb), f"File {fixed_pdb} does not exist"
        # if it is, copy it to the tmp folder we deal with
        shutil.copy(fixed_pdb, f'tmp/{filename}/{filename}_fixed.pdb')

        # do the same copying for _solvated.pdb
        solvated_pdb = os.path.join(restart_rfile_path, f'{filename}_solvated.pdb')
        assert os.path.exists(solvated_pdb), f"File {solvated_pdb} does not exist"
        shutil.copy(solvated_pdb, f'tmp/{filename}/{filename}_solvated.pdb')

        # TODO: do I need to copy the kernels? the colvar? anything like that?

        restart_checkpoint = os.path.join(restart_rfile_path, f'{filename}.chk')
        assert os.path.exists(restart_checkpoint), f"File {restart_checkpoint} does not exist"


        def extract_contact_pairs_str(plumed_file):
            import re
            with open(plumed_file, 'r') as f:
                content = f.read()
            # Find all matches
            matches = re.findall(r'ATOMS\d+=\d+,\d+', content)

            # Join matches into a single string with newlines
            result_string = '\n\t'.join(matches)
            result_string = f"	{result_string}"  # Add leading tab for formatting
            
            return result_string
        
        contact_pairs_str = extract_contact_pairs_str(plumed_file)

    from openmm.unit import nanoseconds, picoseconds
    chk_interval = int((1 * nanoseconds) / (timestep * 0.001 * picoseconds))

    logger.info(f"Checkpoint interval every {chk_interval} steps")

    # create the input for OPES
    from src.plumed.io import create_opes_input
    temperature = 300
    config = {
        'pace': pace,
        'barrier': barrier,
        'temperature': temperature,
        'stride': pace,
        'cutoff': CUTOFF,
        'restart_rfile': restart_rfile,
        'state_wstride': chk_interval,
    }
    create_opes_input(
        filepath=filepath, 
        cv_string=contact_pairs_str,
        config=config,
        type='opes',
        output_dir=output_dir
        )

    from src.plumed.opes import opes
    opes(
        filename=filename, 
        mdtime=mdtime, 
        device_index=str(device_index),
        timestep=timestep,
        temperature=temperature,
        restart_checkpoint=restart_checkpoint,
        device=device,
        output_dir=output_dir,
        chk_interval=chk_interval
        )

def restart(filename=None):
    if filename is None:
        raise ValueError('Filename is required')

    # just run stability from a checkpoint
    print("-----Running stability.py-----")
    from src.stability import stability
    # HACK: set to 90 ns, since first run was 10 ns
    stability(filename=filename, mdtime=90, restart=True)


import fire
if __name__ == '__main__':
    fire.Fire(main)