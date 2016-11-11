# -*- coding: utf-8 -*-
"""
Radation and temperature time series MEGAN

@author: Jordan Capnerhurst 2016
"""

# import things
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from pylab import rcParams

# Set colour map
cool = cm = plt.get_cmap('nipy_spectral')

f = netCDF4.Dataset('O:/Honours_data/JMEGAN/New Jeremy data/rad&t/rad1.nc',
                    'r')


# [ time, source, lon, lat ]
rad = (f.variables['GSW'][:, 14:24, 0, :, :]) # Radiation at ground level
rad2 = (f.variables['GSW'][:, 0:14, 0, :, :]) 
temp = (f.variables['TEMP2'][:, 14:24, 0, :, :]) # temp at 2m
temp2 = (f.variables['TEMP2'][:, 0:14, 0, :, :])

# concatenate back to 24 hour period
vcon = np.concatenate((rad, rad2), axis=1)
tercon = np.concatenate((temp, temp2), axis=1)

# Reshape array for month
rrad = vcon.reshape(672, 128, 153)
rtemp = tercon.reshape(672, 128, 153)-273

#calculate standard devation along lat and lons
stdev1 = np.std(rrad, axis=(1, 2))
stdev2 = np.std(rtemp, axis=(1, 2))

# Averages
v1 = rrad.mean(axis=(1, 2))
ter1 = rtemp.mean(axis=(1, 2))

# X axis dates instead of times
date = np.arange(v1.shape[0])  # assume that delta time between data is 1
date21 = (date/24.)  # use days instead of hours

# change plot size for monthly
rcParams['figure.figsize'] = 15, 5

plt.plot(date21, ter1 , linestyle='-', linewidth=2.0, c='k', label='Ambient Temperature')

# Create standard devation fill
plt.fill_between(date21, ter1-stdev2, ter1+stdev2, alpha=0.7, edgecolor='black', 
                 facecolor='cyan', linewidth=0.3, label='1 Standard devation' )

# asthetics
plt.xlabel('Day')
plt.ylabel('Temperature ($^o$c)')
plt.title('Monthly Average Temperature February 2011 from MEGAN')
#plt.ylim(14, 36)
plt.xlim(0, 27)

# ticks
#plt.xticks(range(0, 28, 1), [str(i) for i in range(0, 28, 1)])
#plt.yticks(np.arange(min(ter1)-5.24, max(ter1 )+11, 2))  # use floats for y 

# Legend
plt.legend(loc='upper right')

# change plot size for daily
rcParams['figure.figsize'] = 7, 4

plt.show()
