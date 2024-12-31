import shutil
import logging

from src.utils import get_checkpoint_interval
from src.analysis.utils import get_file_by_extension
from src.plumed.utils import get_pace_from_metad, prepare_plumed_file_for_restart, read_openmm_output

logger = logging.getLogger(__name__)

def get_last_checkpoint_step(out_file, timestep):
    df_out = read_openmm_output(out_file)

    num_steps_per_checkpoint = get_checkpoint_interval(timestep)
    time_per_checkpoint = num_steps_per_checkpoint * timestep * 0.001 # in ps

    # Go through the dataframe in reverse order to find the first time that's a multiple of time_per_checkpoint
    for idx in reversed(df_out.index):
        current_time = df_out.loc[idx, 'Time']
        if round(current_time, 0) % time_per_checkpoint == 0:
            return int(df_out.loc[idx, 'Step'])
    
    raise ValueError("No corresponding timestep found in the output file for the checkpoint timestep.")

def process_hills_for_restart(
    input_dir: str,
    output_dir: str,
    filename: str,
    timestep: float,
    metad_pace: int
) -> None:
    """Process HILLS file for restart run and write new HILLS file.
    
    Args:
        hills_file: Path to the HILLS file
        out_file: Path to the OpenMM output file
        output_dir: Directory to write the new HILLS file
        filename: Base filename for the new HILLS file
        timestep: Simulation timestep in fs
        metad_pace: MetaD pace in steps
    """
    hills_file = get_file_by_extension(input_dir, '.hills')
    out_file = get_file_by_extension(input_dir, '.out')

    # Calculate the last hill time based on the last checkpoint
    last_checkpoint_step = get_last_checkpoint_step(out_file, timestep)
    num_hills_before_checkpoint = last_checkpoint_step // metad_pace
    last_hill_time = int(num_hills_before_checkpoint * metad_pace * timestep * 0.001)  # in ps
    
    # Process the HILLS file
    with open(hills_file, 'r') as f:
        lines = f.readlines()

    # Get the index of the first line of the last restart
    # Do it by finding the 3 from the end line with #
    counter = 0
    for index, line in reversed(list(enumerate(lines))):
        if line.startswith('#'):
            counter += 1
        if counter == 3:
            break
    # i is the index of the first line of the last restart header; we need to modify this on
    first_index_to_modify = index + 3 # this is the first index of the last restart, header size is 3
    
    # Keep header lines (starting with #) and hills up to time_of_last_hill
    filtered_lines = []
    for index, line in enumerate(lines):
        if index < first_index_to_modify:
            filtered_lines.append(line)
            continue

        # Hills file columns are space-separated, time is in the first column
        try:
            time = float(line.split()[0])
            if time <= last_hill_time:
                filtered_lines.append(line)
        except (ValueError, IndexError):
            continue

    # Write the filtered lines to the new HILLS file
    with open(f'{output_dir}/{filename}.hills', 'w') as f:
        f.writelines(filtered_lines)
    
    # Also, re-write the original file, with the cut-offed post-checkpoint hills
    with open(hills_file, 'w') as f:
        f.writelines(filtered_lines)

    # Also cut the colvar accordingly
    from src.analysis.io import read_colvar_file
    colvar_file = get_file_by_extension(input_dir, '.colvar')
    colvar_df = read_colvar_file(colvar_file)
    colvar_df = colvar_df[colvar_df['time'] <= last_hill_time]
    indices_to_keep = 3 + len(colvar_df) # 3 is the header size
    with open(colvar_file, 'r') as f:
        lines = f.readlines()
    with open(colvar_file, 'w') as f:
        f.writelines(lines[:indices_to_keep])

# move hills file to output dir
# chop off lines that are after the checkpoint
# get the checkpoint file
# figure out the timestep
# HILLS file numbers every deposited kernel
# .out prints every 100 ps (i.e. 50,000 time steps)
# we finished at (were WALLtimed) at 221,900 ps, i.e. 221.9 ns
# how often do we get a checkpoint?
# we save every 1 nanoseconds, according to get_checkpoint_interval
def validate_metad_pace(input_dir: str, config: dict) -> None:
    """Validate that the MetaD pace matches between previous and current run."""
    plumed_file = get_file_by_extension(input_dir, 'plumed.dat')
    metad_pace = get_pace_from_metad(plumed_file)
    assert metad_pace == config['metad.pace'], 'Previous MetaD pace does not match the current one.'

def prepare_restart_files(input_dir: str, output_dir: str, filename: str) -> None:
    """Prepare PLUMED input file for restart."""
    plumed_file = get_file_by_extension(input_dir, 'plumed.dat')
    prepare_plumed_file_for_restart(plumed_file, output_dir, filename)

def setup_metad_restart(
    input_dir: str,
    output_dir: str,
    filename: str,
    timestep: float,
    config: dict
) -> str:
    """Main function to set up MetaD restart.
    
    Returns:
        str: Path to the restart checkpoint file
    """
    try:
        validate_metad_pace(input_dir, config)
        prepare_restart_files(input_dir, output_dir, filename)
        process_hills_for_restart(input_dir, output_dir, filename, timestep, config['metad.pace'])
        
        logger.info("Restarting MetaD as requested...")
        return get_file_by_extension(input_dir, '.chk')
    except Exception as e:
        logger.error(f"Error during MetaD restart setup: {e}")
        shutil.rmtree(output_dir)
        raise e


# TODO: this is a mess, needs to be cleaned up


# print("Getting CVs from previous plumed.dat")
# # assume it's in the same folder as restart_rfile
# restart_rfile_path = os.path.dirname(restart_rfile)
# plumed_file = os.path.join(restart_rfile_path, f'{filename}_plumed.dat')
# assert os.path.exists(plumed_file), f"File {plumed_file} does not exist"
# shutil.copy(plumed_file, f"{output_dir}/{filename}_plumed.dat")

# # assume fixed.pdb also in the same folder
# fixed_pdb = os.path.join(restart_rfile_path, f'{filename}_fixed.pdb')
# assert os.path.exists(fixed_pdb), f"File {fixed_pdb} does not exist"
# # if it is, copy it to the tmp folder we deal with
# shutil.copy(fixed_pdb, f'{output_dir}/{filename}_fixed.pdb')

# # do the same copying for _solvated.pdb
# solvated_pdb = os.path.join(restart_rfile_path, f'{filename}_solvated.pdb')
# assert os.path.exists(solvated_pdb), f"File {solvated_pdb} does not exist"
# shutil.copy(solvated_pdb, f'{output_dir}/{filename}_solvated.pdb')

# # TODO: do I need to copy the kernels? the colvar? anything like that?

# restart_checkpoint = os.path.join(restart_rfile_path, f'{filename}.chk')
# assert os.path.exists(restart_checkpoint), f"File {restart_checkpoint} does not exist"

# def extract_contact_pairs_str(plumed_file):
#     import re
#     with open(plumed_file, 'r') as f:
#         content = f.read()
#     # Find all matches
#     matches = re.findall(r'ATOMS\d+=\d+,\d+', content)

#     # Join matches into a single string with newlines
#     result_string = '\n\t'.join(matches)
#     result_string = f"	{result_string}"  # Add leading tab for formatting
    
#     return result_string

# contact_pairs_str = extract_contact_pairs_str(plumed_file)
