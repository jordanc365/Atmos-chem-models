# -*- coding: utf-8 -*-
"""
Create emissions Timeseries MEGAN WRF

@author: Jordan Capnerhurst 2016
"""

# import things
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from pylab import rcParams

# Set colour map
cool = cm = plt.get_cmap('nipy_spectral')

f = netCDF4.Dataset('O:/Honours_data/JMEGAN/FebMegD3.nc', 'r')


# [ time, source, lon, lat ]
v = (f.variables['ISOP'][1:29, 15:24, 0, :, 10:125]*68.12)
v2 = (f.variables['ISOP'][1:29, 0:15, 0, :, 10:125]*68.12)
ter = (f.variables['TERP'][1:29, 15:24, 0, :, 10:125]*136.298)
ter2 = (f.variables['TERP'][1:29, 0:15, 0, :, 10:125]*136.298)

# concatenate back to 24 hour period
vcon = np.concatenate((v, v2), axis=1)
tercon = np.concatenate((ter, ter2), axis=1)

# Add emissions along all arrays
total = (np.sum(vcon, axis=(2, 3))*3.6)
total1 = (np.sum(tercon, axis=(2, 3))*3.6)

# Averages
v1 = total.mean(axis=(0))
ter1 = total1.mean(axis=(0))

# standardss devitations
stdevjmeg1 = np.std(total, axis=(0))
stdevjmeg2 = np.std(total1, axis=(0))



# Add monoterpenes and isoprene to same array daily
totd = [sum(x) for x in zip(stdevjmeg1, stdevjmeg2)]
stdtot = np.array(totd)  # milimoles to kilograms

# Standard deviations
totd = [sum(x) for x in zip(v1, ter1)]
totarrd = np.array(totd)  # milimoles to kilograms


# X axis dates instead of times
#date = np.arange(rv1.shape[0])  # assume that delta time between data is 1
#date21 = (date/24.)  # use days instead of hours

# change plot size for monthly
rcParams['figure.figsize'] = 20, 5
#plt.plot(date21, rter1, linestyle='-', linewidth=3.0, c='g', label=' MEGAN  Emissions')
plt.plot(totarrd, linestyle='-', linewidth=3.0, c='b', label=' CTM Biogenic Emissions')
 #KOEH Data
# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temperature')
#plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temperature')

# asthetics
plt.xlabel('Day')
plt.ylabel('Isoprene Biogenic Emissions (kg)')
plt.title('Isoprene Emissions March 2011')
plt.ylim(0, 25)
plt.xlim(0, 27)


# ticks
plt.xticks(range(1, 28, 1), [str(i) for i in range(1, 28, 1)])
plt.yticks(np.arange(min(totarr2)-541, max(totarr2)+50000, 50000))  # use floats for y 

# Legend
# plt.legend(['MEGAN', 'MEGAN Coupled to CTM', 'Skin Temp'], loc='upper right')

# change plot size for daily
rcParams['figure.figsize'] = 7, 4

plt.show()
