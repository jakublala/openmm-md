from openmmtools.multistate import ParallelTemperingAnalyzer
import pathlib

def main():
    storage_path = pathlib.Path('replica_exchange.nc')
    import netCDF4
    nc = netCDF4.Dataset(storage_path, 'r')
    print(nc.variables.keys())
    n_replicas = nc['states'].shape[0]
    n_iterations = nc['states'].shape[1]

    from openmm.app import PDBFile
    from openmm.unit import nanometer, angstrom
    pdb_path = '../../data/241010_FoldingUponBinding/output/SUMO-1C/241128-MetaD/sumo1c_equilibrated.pdb'
    pdb = PDBFile(pdb_path)



    # loading traj and write to dcd format
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader
    coord_traj = np.zeros((n_replicas, 20, pdb.topology.getNumAtoms(), 3)) * nanometer
    for time in range(0, n_iterations, 10):
        print(nc['metadata'])
        # state_samplers = reporter.read_sampler_states(time)
        for i, state_sampler in enumerate(state_samplers):
            coord_traj[i, int(time/10)] += state_sampler.positions

    # only taking the first replica trajectory
    u = mda.Universe(pdb_path, coord_traj[0].value_in_unit(angstrom), format=MemoryReader)
    with mda.Writer("test.dcd", u.select_atoms("protein").n_atoms) as W:
        for ts in u.trajectory:
            W.write(u.select_atoms("protein"))


if __name__ == "__main__":
    main()