


if __name__ == '__main__':
    import netCDF4
    import numpy as np
    import matplotlib.pyplot as plt
    import pdb

    import os


    nc = netCDF4.Dataset('tmp/replica_exchange.nc', 'r')

    print(nc.variables.keys())
    print(nc.variables['states'][:])


    # get energies
    energies = nc.variables['energies'][:]
    
    # print(energies[:, 0, :])
    # plot these energies, where x-axis is energy, and y-axis is population density, histogram
    plt.hist(energies[:, 0, :], bins=100)
    plt.savefig('tmp/energies.png')
