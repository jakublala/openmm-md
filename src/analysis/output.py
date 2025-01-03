
import pandas as pd
from src.utils import get_restarted_files_by_extension
import os 

# HACK
LOGGING_INTERVAL = 100 # ps per frame
CHECKPOINT_INTERVAL = 1000 # ps per checkpoint, 1 ns per checkpoint

def stitch_openmm_outputs(directory, system, date):
    try:
        os.remove(f"{directory}/{system}_full.out")
    except FileNotFoundError:
        pass

    files = get_restarted_files_by_extension(directory, '.out')
    
    assert files == sorted(files), "Files are not sorted"

    df = pd.DataFrame()
    for i, file in enumerate(files):
        _df = pd.read_csv(file, sep='\t')
        # Get last index that's a multiple of 10 (checkpoint occurs every 10 frames)
        # don't do this for the last file
        if i < len(files) - 1:
            last_checkpoint_idx = len(_df) - (len(_df) % (CHECKPOINT_INTERVAL // LOGGING_INTERVAL))
            _df = _df.iloc[:last_checkpoint_idx]
        df = pd.concat([df, _df])

    time_col = df['Time (ps)'].to_list()
    time_col = [round(t, 0) for t in time_col]
    assert all(t % LOGGING_INTERVAL == 0 for t in time_col), "Time column is not a multiple of LOGGING_INTERVAL"
    assert all(time_col[i+1] - time_col[i] == LOGGING_INTERVAL for i in range(len(time_col) - 1)), "Consecutive entries are not LOGGING_INTERVAL apart"

    df.to_csv(f"{directory}/{system}_full.out", sep='\t', index=False)
