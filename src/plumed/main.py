import sys
print(sys.executable)
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
import os
import subprocess

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
        timestep=4
        ):
    if filepath is None:
        raise ValueError('Filepath is required')

    filename = os.path.basename(filepath).split('.')[0]

    print(f'==================== Running {filename} ====================')
    
    if not os.path.exists('tmp'):
        # create the tmp folder
        run_command('mkdir tmp')


    # 1. load the PDB and fix errors
    from ..fixer import fixer
    fixer(filepath=filepath)

    # 1.5 get OPES preparation
    from .cv import get_interface_contact_indices
    contact_indices = get_interface_contact_indices(filepath)
    print(contact_indices)

    assert 0 == 1, "Stop here"
    # 2. minimize the structure with LBFGS and H atoms mobile
    from ..relax import minimize
    minimize(
        filename=filename, 
        max_iterations=0, 
        device_index=str(device_index),
        constraints=None
        )
    
    try:
        print("-----Running stability.py-----")
        now = datetime.now()
        dt_string = now.strftime("%y%m%d_%H%M%S")

        from src.stability import stability
        stability(
            filename=filename, 
            mdtime=mdtime, 
            device_index=str(device_index),
            timestep=timestep
            )
        
    except Exception as e:
        print(f"Error running stability: {e}")
        if os.path.exists(f'tmp/{filename}.xyz'):
            run_command(f'mv tmp/{filename}.xyz output/{filename}_{dt_string}.xyz')
        if os.path.exists(f'tmp/{filename}.out'):
            run_command(f'mv tmp/{filename}.out output/{filename}_{dt_string}.out')
        if os.path.exists(f'tmp/{filename}.chk'):
            run_command(f'mv tmp/{filename}.chk output/{filename}_{dt_string}.chk')
        if os.path.exists(f"tmp/{filename}_solvated.pdb"):
            run_command(f"mv tmp/{filename}_solvated.pdb output/{filename}_{dt_string}_solvated.pdb")

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