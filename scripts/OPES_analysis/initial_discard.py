from src.analysis.fes import compute_fes
import os
import logging
from joblib import Parallel, delayed


logger = logging.getLogger(__name__)

from src.analysis.colvar import read_colvar_file
from src.analysis.kernels import get_sigmas
from src.analysis.plot import plot_2d_fes
            
# Process each discard range in parallel
from concurrent.futures import ThreadPoolExecutor

            
def process_discard(discard, colvar_df, sigmas):
    logger.info(f"Processing discard {discard} of the trajectory.")
    # Create copy of dataframe with discarded steps
    num_steps = len(colvar_df)
    num_steps_to_discard = int(num_steps * discard) if discard is not None else 0
    df_subset = colvar_df.iloc[num_steps_to_discard:].copy()
    
    # Compute FES for this subset
    cv1_bins, cv2_bins, fes = compute_fes(
        df_subset,
        sigma=sigmas,
        temp=300,
        cvs=['cmap', 'd'],
        outfile=None,
        bias=['opes.bias', 'uwall.bias'],
        n_bins=200
    )
    return cv1_bins, cv2_bins, fes

from src.analysis.plot import plot_all_fes_from_data
def run(date, system, num_runs):
    # discard first 1/10, 1/5, 1/4, 1/3, and 1/2 of the trajectory
    discard_ranges = [None, 1/10, 1/5, 1/4, 1/3]

    cvs = ['cmap', 'd']
    target, binder = system.split("_")
    for run in range(1, num_runs + 1):
        results = []
        # First check if directory exists
        directory = f"../../data/241010_FoldingUponBinding/output/{date}/{target}/{binder}_{run}"
        if not os.path.exists(directory):
            logger.warning(f"System {system} does not exist for experiment {date}")
            continue

        if target == "SUMO":
            target = "sumo"
        
        # Read colvar file once
        colvar_df = read_colvar_file(f"{directory}/{target}_{binder}.colvar")
        
        # Get sigmas once
        sigmas = get_sigmas(directory, target, binder, run, cvs)
        
        logger.info(f"Processing {system} run {run}")
        for discard in discard_ranges:
            result_fes = process_discard(
                discard, 
                colvar_df, 
                sigmas
            )
            results.append(result_fes)

        # make a subplot next to one another of num_runs
        outfile = f"{directory}/{system}_initial_discard.png"
        plot_all_fes_from_data(results, outfile, target, binder, cvs, labels=[f"{discard=}" for discard in discard_ranges])
            




def main(system):
    num_runs = 5
    # systems = [
    #     'A-synuclein_alpha', 'A-synuclein_general', 
    #     'CD28_alpha', 'CD28_beta', 'CD28_partial', 'CD28_general',
    #     'p53_1', 'p53_2', 'p53_end',
    #     'sumo_1', 'sumo_1c'
    # ]
    date = "241029"
    run(date, system, num_runs)


import fire
if __name__ == "__main__":
    fire.Fire(main)
