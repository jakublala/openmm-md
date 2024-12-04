import glob
def get_file_by_extension(directory, extension, assert_exists=True):
    """Find first file with given extension in directory.
    
    Args:
        directory: Directory to search in
        extension: File extension to look for (e.g. '.colvar', '.h5')
        
    Returns:
        str: Path to first matching file
        
    Raises:
        FileNotFoundError: If no file with extension is found
    """
    files = glob.glob(f"{directory}/*{extension}")
    for i, file in enumerate(files):
        if 'bck.' in file:
            logger.warning(f"Found a backup file ({file}), ignoring it...")
            files.pop(i)
    
    if len(files) == 0:
        if assert_exists:
            raise FileNotFoundError(f"No file with extension {extension} found in {directory}")
        else:
            return False
    return files[0]

from src.constants import kB
def convert_deltaG_to_kBT(deltaG_kJmol, TEMP):
    return deltaG_kJmol / (kB * TEMP)


import os
import logging
logger = logging.getLogger(__name__)
def fix_fucked_up_naming(project, system, date):
    if '-' in str(date):
        # check if this folder or the one with _ exists and adjust accordingly
        if os.path.exists(f"../../data/{project}/output/{system}/{date}"):
            pass
        elif os.path.exists(f"../../data/{project}/output/{system}/{date.replace('-', '_')}"):    
            date = date.replace('-', '_')
            logger.warning(f"Replacing '-' with '_' in date: {date}")
            logger.warning("Future simulations should avoid '-' in date")
        else:
            raise ValueError(f"No directory found for {date} or {date.replace('-', '_')}")
    return date
