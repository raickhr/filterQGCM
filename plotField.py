from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

fileName = 'recNo_0001.filterLen_0100.nc'
fieldName = 'EddyPowerPerArea'
ds = Dataset(fileName)

ff = np.array(ds.variables[fieldName])

plt.pcolormesh(ff[0,:,:])
plt.colorbar()
plt.show()
