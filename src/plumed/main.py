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
        timestep=2
        ):
    if filepath is None:
        raise ValueError('Filepath is required')

    filename = os.path.basename(filepath).split('.')[0]

    print(f'==================== Running {filename} ====================')
    os.makedirs(f'tmp/{filename}', exist_ok=True)


    # 1. load the PDB and fix errors
    from src.fixer import fixer
    fixer(filepath=filepath)


    # 2. minimize the structure with LBFGS and H atoms mobile
    from src.relax import minimize
    minimize(
        filename=filename, 
        max_iterations=0, 
        device_index=str(device_index),
        constraints=None
        )
    
    # 2.5 get OPES preparation
    from src.plumed.cv import get_interface_contact_indices
    cutoff = 0.8
    contact_indices = get_interface_contact_indices(filename=filename, cutoff=cutoff)

    contact_pairs_str = ""
    for i, pair in enumerate(contact_indices):
        if i != 0:
            contact_pairs_str += f"\n\tATOMS{i+1}={pair[0]},{pair[1]}"
        else:
            contact_pairs_str += f"\tATOMS{i+1}={pair[0]},{pair[1]}"


    # create the input for OPES
    from src.plumed.io import create_opes_input
    temperature = 300
    config = {
        'pace': 500,
        'barrier': 50,
        'temperature': temperature,
        'stride': 500,
        'cutoff': cutoff,
        'restart_rfile': restart_rfile,
    }
    create_opes_input(
        filepath=filepath, 
        cv_string=contact_pairs_str,
        config=config
        )

    now = datetime.now()
    dt_string = now.strftime("%y%m%d_%H%M%S")

    from src.plumed.opes import opes
    opes(
        filename=filename, 
        mdtime=mdtime, 
        device_index=str(device_index),
        timestep=timestep,
        temperature=temperature
        )
        
    # except Exception as e:
    #     print(f"Error running stability: {e}")
    #     if os.path.exists(f'tmp/{filename}.xyz'):
    #         run_command(f'mv tmp/{filename}.xyz output/{filename}_{dt_string}.xyz')
    #     if os.path.exists(f'tmp/{filename}.out'):
    #         run_command(f'mv tmp/{filename}.out output/{filename}_{dt_string}.out')
    #     if os.path.exists(f'tmp/{filename}.chk'):
    #         run_command(f'mv tmp/{filename}.chk output/{filename}_{dt_string}.chk')
    #     if os.path.exists(f"tmp/{filename}_solvated.pdb"):
    #         run_command(f"mv tmp/{filename}_solvated.pdb output/{filename}_{dt_string}_solvated.pdb")

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