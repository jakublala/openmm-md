import sys
print(sys.executable)

import subprocess

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        print(f"Error running command: {command}")
        print(f"Return code: {process.returncode}")
        print(f"Standard output:\n{stdout.decode()}")
        print(f"Standard error:\n{stderr.decode()}")
        sys.exit(1)  # Exit the script with an error code
    
    return stdout.decode()

def main(filename=None):
    if filename is None:
        raise ValueError('Filename is required')
    
    print("-----Running fixer.py-----")
    output = run_command(f'python fixer.py --filename {filename}')
    print(output)
    
    print("-----Running minimize.py-----")
    output = run_command(f'python minimize.py --filename {filename}')
    print(output)

    print("-----Running stability.py-----")
    output = run_command(f'python stability.py --filename {filename} --nsteps=5000')
    print(output)


if __name__ == '__main__':
    # filenames = ['S1_Best_A', 'S1_Best_AB', 'S2_Best_A', 'S2_Best_AB']
    filenames = ['tutorial_example']
    for filename in filenames:
        print(f'==================== Running {filename} ====================')
        main(filename=filename)
    
    print("All scripts completed successfully.")