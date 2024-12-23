


if __name__ == '__main__':
    import netCDF4
    import numpy as np
    import matplotlib.pyplot as plt
    import pdb

    import os


    nc = netCDF4.Dataset('tmp/replica_exchange.nc', 'r')

    print(nc.variables.keys())
    print(nc.variables['states'][:])
