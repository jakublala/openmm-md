from openmm.unit import nanoseconds, picoseconds
def get_checkpoint_interval(timestep):
    return int((1 * nanoseconds) / (timestep * 0.001 * picoseconds)) 