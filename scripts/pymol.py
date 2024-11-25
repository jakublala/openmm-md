from pymol import cmd

def translate_trajectory(object_name, dx=None, dy=None, dz=None):
    """
    Translates an object across all frames in the trajectory, automatically 
    retrieving unit cell dimensions if dx, dy, or dz is not provided.

    Parameters:
    - object_name: Name of the object to translate
    - dx, dy, dz: Translation values along the x, y, and z axes. If None, these
                  are computed as half the unit cell dimensions.
    """
    # Get symmetry info
    symmetry_info = cmd.get_symmetry(object_name)
    if symmetry_info is None:
        raise ValueError(f"Object '{object_name}' does not have symmetry information.")

    # Extract unit cell dimensions
    a, b, c = symmetry_info[0], symmetry_info[1], symmetry_info[2]

    # Default translations: move by half the unit cell dimensions
    dx = dx if dx is not None else a / 2
    dy = dy if dy is not None else b / 2
    dz = dz if dz is not None else c / 2

    print(f"Translating '{object_name}' by: dx={dx}, dy={dy}, dz={dz}")

    # Apply translation to all frames
    num_states = cmd.count_states(object_name)
    for state in range(1, num_states + 1):  # PyMOL states are 1-indexed
        cmd.alter_state(state, object_name,
                        f"(x, y, z) = (x + {dx}, y + {dy}, z + {dz})")
    cmd.refresh()
    print(f"Translation complete for {num_states} frames of '{object_name}'.")