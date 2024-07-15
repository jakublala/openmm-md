import sys
print(sys.executable)

import os

def main(
        filename=None
):
    if filename is None:
        raise ValueError('Filename is required')
    
    print("-----Running fixer.py-----")
    os.system(f'python fixer.py --filename {filename}')
    
    print("-----Running minimize.py-----")
    os.system(f'python minimize.py --filename {filename}')

    print("-----Running stability.py-----")
    os.system(f'python stability.py --filename {filename}')


if __name__ == '__main__':
    # filenames = ['S1_Best_A', 'S1_Best_AB', 'S2_Best_A', 'S2_Best_AB']
    filenames = ['S1_Best_AB', 'S2_Best_A', 'S2_Best_AB']
    for filename in filenames:
        print(f'==================== Running {filename} ====================')
        main(filename=filename)