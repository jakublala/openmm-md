import sys
print(sys.executable)

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

def main(filename=None, device_index=0):
    if filename is None:
        raise ValueError('Filename is required')
    
    print(f'==================== Running {filename} ====================')
    
    if not os.path.exists('tmp'):
        # create the tmp folder
        run_command('mkdir tmp')


    # print("-----Running fixer.py-----")
    from fixer import fixer
    fixer(filename=filename)


    # do 10 attempts
    for i in range(10):
        print("-----Running minimize.py-----")
        from minimize import minimize
        minimize(filename=filename, max_iterations=1000, device_index=str(device_index))

        try:
            print("-----Running stability.py-----")
            print(f"Running for the {i+1}th time...")
            now = datetime.now()
            dt_string = now.strftime("%y%m%d_%H%M%S")

            from stability import stability
            stability(filename=filename, mdtime=100, device_index=str(device_index))
            break
        except Exception as e:
            print(f"Error running stability: {e}")
            print("Trying again...")
            # remove the xtc file
            # if this file exists run the command
            if os.path.exists(f'tmp/{filename}.xtc'):
                run_command(f'mv tmp/{filename}.xtc output/{filename}_{dt_string}.xtc')
            if os.path.exists(f'tmp/{filename}.out'):
                run_command(f'mv tmp/{filename}.out output/{filename}_{dt_string}.out')
            if os.path.exists(f'tmp/{filename}.chk'):
                run_command(f'mv tmp/{filename}.chk output/{filename}_{dt_string}.chk')
            if os.path.exists(f"tmp/{filename}_solvated.pdb"):
                run_command(f"mv tmp/{filename}_solvated.pdb output/{filename}_{dt_string}_solvated.pdb")
            continue

def restart(filename=None):
    if filename is None:
        raise ValueError('Filename is required')

    # just run stability from a checkpoint
    print("-----Running stability.py-----")
    from stability import stability
    # HACK: set to 90 ns, since first run was 10 ns
    stability(filename=filename, mdtime=90, restart=True)


import fire
if __name__ == '__main__':
    fire.Fire(main)