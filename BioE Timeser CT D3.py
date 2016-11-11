# -*- coding: utf-8 -*-
"""
Create Timeseries CTM 2011

@author: Jordan Capnerhurst 2016
"""

# import netCDF
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from pylab import rcParams
fc = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc', 'r')

# plot Daily average
vbios = fc.variables['store_Bio'][:, 14:24, 0, 6:46, 6:48]
vbios2 = fc.variables['store_Bio'][:, 0:14, 0, 6:46, 6:48]
t = fc.variables['temp_a'][:, 14:24, 0, :, :]
t2 = fc.variables['temp_a'][:, 0:14, 0, :, :]
g = fc.variables['skin_temp'][:, :, :, :]

# concatenate back to 24 hour period
vcon = np.concatenate((vbios, vbios2), axis=1)
tempcon = np.concatenate((t, t2), axis=1)

totalz = (np.sum(vcon, axis=(2, 3)))*0.91

# Mean over array
bios1 = vcon.mean(axis=(2, 3))/9
t1 = tempcon.mean(axis=(2, 3))
g1 = g.mean(axis=(2, 3))

# Reshape array for month
biosr2 = np.reshape(t1, 648)  # Due to Dimensions being 3x3
#r3 = np.reshape(t1, 648)
#r4 = np.reshape(g1, 675)


# X axis dates instead of times for monthly 
date = np.arange(biosr2.shape[0])  # assume  delta time between data is 1 hour
date1 = (date/24.)  # use days instead of hours
ar = np.array(date1)

# change plot size for monthly
rcParams['figure.figsize'] = 15, 5

plt.plot(ar, biosr2, linestyle='--', linewidth=3.0, c='b',
         label='Biogenic Emissions')


# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temp')
# plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temp')

# asthetics
plt.xlabel('Day')
plt.ylabel('Temperature ($^o$c)')
plt.title('Monthly Ambiant temperature Feburary 2011 from CTM')
plt.ylim(14, 35)
plt.xlim(0, 27)


# tick
plt.xticks(range(0, 28, 1), [str(i) for i in range(0, 28, 1)])
plt.yticks(np.arange(min(biosr2)-1.45, max(biosr2)+2, 1))  # use floats for y value

# Legend
#plt.legend(['Store Bio', 'Ambiant temp', 'Skin Temp'], loc='upper left')

# change plot size
rcParams['figure.figsize'] = 7, 4

plt.show()
# plt.savefig('./bmap_syd.png')
