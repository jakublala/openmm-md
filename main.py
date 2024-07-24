import sys
print(sys.executable)

import os
import subprocess

# current datetime
from datetime import datetime
now = datetime.now()
dt_string = now.strftime("%y%m%d_%H%M%S")

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

def main(filename=None):
    if filename is None:
        raise ValueError('Filename is required')
    
    # if tmp folder exists, remove it
    if os.path.exists('tmp'):
        run_command('rm -r tmp')
    if not os.path.exists('tmp'):
        # create the tmp folder
        run_command('mkdir tmp')


    # print("-----Running fixer.py-----")
    from fixer import fixer
    fixer(filename=filename)


    # do 10 attempts
    for i in range(100):
        print("-----Running minimize.py-----")
        from minimize import minimize
        minimize(filename=filename, max_iterations=1000)

        try:
            print("-----Running stability.py-----")
            print(f"Running for the {i+1}th time...")
            from stability import stability
            stability(filename=filename, mdtime=10)
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
            continue

def restart(filename=None):
    if filename is None:
        raise ValueError('Filename is required')

    # just run stability from a checkpoint
    print("-----Running stability.py-----")
    from stability import stability
    # HACK: set to 90 ns, since first run was 10 ns
    stability(filename=filename, mdtime=90, restart=True)


if __name__ == '__main__':
    # filenames = ['S1_Best_A', 'S1_Best_AB', 'S2_Best_A', 'S2_Best_AB']
    filenames = ['S1_Best_AB']
    for filename in filenames:
        print(f'==================== Running {filename} ====================')
        main(filename=filename)
    
    print("All scripts completed successfully.")