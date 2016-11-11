# -*- coding: utf-8 -*-
"""
Create monthly/ daily emission timeseries 3x3 JMEGAN

@author: Jordan Capnerhurst 2016
"""

# import things
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from pylab import rcParams

# Set colour map
cool = cm = plt.get_cmap('nipy_spectral')

f = netCDF4.Dataset('O:/Honours_data/JMEGAN/FebJMEGv2.nc', 'r')


# [ time, source, lon, lat ]
v = (f.variables['ISOP'][1:29, 15:24, 3:46, 6:50]*68.12)
v2 = (f.variables['ISOP'][1:29, 0:15, 3:46, 6:50]*68.12)
ter = (f.variables['TERP'][1:29, 15:24, 3:46, 6:50]*136.298)
ter2 = (f.variables['TERP'][1:29, 0:15, 3:46, 6:50]*136.298)

# concatenate back to 24 hour period
vcon = np.concatenate((v, v2), axis=1)
tercon = np.concatenate((ter, ter2), axis=1)

# Add emissions along all arrays
total = (np.sum(vcon, axis=(2, 3)))
total1 = (np.sum(tercon, axis=(2, 3)))

# Averages
v1 = vcon.mean(axis=(0, 2, 3))
ter1 = tercon.mean(axis=(0, 2, 3))


# Reshape array for month
rv1 = np.reshape(v1, 648)
rter1 = np.reshape(ter1, 648)


# Add monoterpenes and isoprene to same array Monthly
tot2 = [sum(x) for x in zip(rv1, rter1)]
totarr2 = np.array(tot2)*3.6  # milimoles to kilograms

# Add monoterpenes and isoprene to same array daily
totd = [sum(x) for x in zip(v1, ter1)]
totarrd = np.array(totd)*3.6/9# milimoles to kilograms


# X axis dates instead of times
date = np.arange(totarr2.shape[0])  # assume that delta time between data is 1
date21 = (date/24.)  # use days instead of hours

# change plot size for monthly
rcParams['figure.figsize'] = 25, 5

plt.plot(date21, totarr2, linestyle='-', linewidth=3.0, c='g', label=' MEGAN  Emissions')
#plt.plot(v1, linestyle='-', linewidth=3.0, c='b', label=' CTM Biogenic Emissions') #KOEH Data
# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temperature')
#plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temperature')

# asthetics
plt.xlabel('Day')
plt.ylabel('Total Biogenic Emissions (kg/$m^2$/hour)')
plt.title('Monoterpene Emissions March 2011')
#plt.ylim(0, 7)
plt.xlim(0, 27)


# ticks
plt.xticks(range(1, 28, 1), [str(i) for i in range(1, 28, 1)])
plt.yticks(np.arange(min(totarr2)-0.067, max(totarr2)+1, 1))  # use floats for y 

# Legend
# plt.legend(['MEGAN', 'MEGAN Coupled to CTM', 'Skin Temp'], loc='upper right')

# change plot size for daily
rcParams['figure.figsize'] = 7, 4

plt.show()
