import os
from tqdm import tqdm

def main(filepath='../../data/input/241010_FUB/A-synuclein/A-synuclein_alpha.pdb'):    
    filename = os.path.basename(filepath).split('.')[0]
    os.makedirs(f'tmp/{filename}', exist_ok=True)
    device_index = 0

    minimized_filepath = f'tmp/{filename}/{filename}_solvated.pdb'

    if not os.path.exists(minimized_filepath):
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
    import numpy as np
    import matplotlib.pyplot as plt
    from src.plumed.cv import get_interface_contact_indices

    cutoffs = np.linspace(0, 3.0, 31)  # 31 points from 0 to 3.0
    contact_counts = []

    for cutoff in tqdm(cutoffs, desc="Calculating contact counts"):
        contact_indices = get_interface_contact_indices(filename=filename, cutoff=cutoff)
        contact_counts.append(len(contact_indices))

    plt.figure(figsize=(10, 6))
    plt.plot(cutoffs, contact_counts, marker='o')
    plt.xlabel('Cutoff (Angstroms)')
    plt.ylabel('Number of Contacts')
    plt.title(f'system = {filename}')
    
    os.makedirs('output', exist_ok=True)
    plt.savefig(f'output/{filename}_contact_cutoff_plot.png', dpi=300)
    plt.close()


def plot_rational_switching_function():
    import matplotlib.pyplot as plt
    import numpy as np
    def rational_function(r, r0, d0, n, m):
        numerator = 1 - ((r - d0) / r0) ** n
        denominator = 1 - ((r - d0) / r0) ** m
        return numerator / denominator

    # Parameters
    r0 = 0.8  # characteristic distance
    d0 = 0.0  # offset
    n = 2     # numerator exponent
    m = 4    # denominator exponent (2n)

    # Generate r values
    r = np.linspace(0, 10, 1000)

    # Calculate function values
    s = rational_function(r, r0, d0, n, m)

    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(r, s)
    plt.grid(True)
    plt.xlabel('r')
    plt.ylabel('s(r)')
    plt.title(f'Rational Switching Function\n r0 = {r0}, d0 = {d0}, n = {n}, m = {m}')
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    plt.savefig('output/rational_switching_function.png', dpi=300)

    plt.show()



    


if __name__ == '__main__':
    plot_rational_switching_function()
    # run this for all systems
    filepaths = [
        '../../data/input/241010_FUB/A-synuclein/A-synuclein_alpha.pdb',
        '../../data/input/241010_FUB/A-synuclein/A-synuclein_general.pdb',
        '../../data/input/241010_FUB/CD28/CD28_alpha.pdb',
        '../../data/input/241010_FUB/CD28/CD28_beta.pdb',
        '../../data/input/241010_FUB/CD28/CD28_general.pdb',
        '../../data/input/241010_FUB/CD28/CD28_partial.pdb',
        '../../data/input/241010_FUB/p53/p53_1.pdb',
        '../../data/input/241010_FUB/p53/p53_2.pdb',
        '../../data/input/241010_FUB/p53/p53_end.pdb',
        '../../data/input/241010_FUB/SUMO/sumo1.pdb',
        '../../data/input/241010_FUB/SUMO/sumo1c.pdb'
    ]

    for filepath in tqdm(filepaths, desc="Processing systems"):
        main(filepath=filepath)