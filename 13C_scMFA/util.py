from pathlib import Path

def construct_path(file_path):
    """
    Checks if a path exists.
    
    Arguments
    ---------
    file_path : str
        path to check
    
    Returns
    -------
    pathlib.Path to the input path.
    """
    path = Path(file_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist.")
    return path