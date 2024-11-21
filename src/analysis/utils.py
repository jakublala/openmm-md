import glob
def get_file_by_extension(directory, extension, assert_exists=True):
    """Find first file with given extension in directory.
    
    Args:
        directory: Directory to search in
        extension: File extension to look for (e.g. '.colvar', '.h5py')
        
    Returns:
        str: Path to first matching file
        
    Raises:
        FileNotFoundError: If no file with extension is found
    """
    files = glob.glob(f"{directory}/*{extension}")
    if len(files) == 0:
        if assert_exists:
            raise FileNotFoundError(f"No file with extension {extension} found in {directory}")
        else:
            return False
    return files[0]