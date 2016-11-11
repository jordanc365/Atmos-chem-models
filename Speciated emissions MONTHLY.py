# -*- coding: utf-8 -*-
"""
Create monthly emissions Timeseries MEGAN WRF

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

# Reshape array for month
rv1 = vcon.reshape(672, 128, 115)
rter1 = tercon.reshape(672,  128, 115 )

# standardss devitations
stdevjmeg1 = np.std(rv1, axis=(1, 2))*3.6
stdevjmeg2 = np.std(rter1, axis=(1, 2))*3.6

# Add emissions along all arrays
#total = (np.sum(vcon, axis=(2, 3))*3.6)
#total1 = (np.sum(tercon, axis=(2, 3))*3.6)

# Averages
v1 = rv1.mean(axis=(1, 2))*3.6
ter1 = rter1.mean(axis=(1, 2))*3.6

# Add monoterpenes and isoprene to same array Monthly
tot2 = [sum(x) for x in zip(v1, ter1)]
totarr2 = np.array(tot2)  # milimoles to kilograms

# stds for total emissions
tot2 = [sum(x) for x in zip(stdevjmeg1, stdevjmeg2)]
totstdjmeg = np.array(tot2)  

# Add monoterpenes and isoprene to same array daily
#totd = [sum(x) for x in zip(v1, ter1)]
#totarrd = np.array(totd)*3.6  # milimoles to kilograms

# X axis dates instead of times
date = np.arange(v1.shape[0])  # assume that delta time between data is 1
date21 = (date/24.)  # use days instead of hours

# change plot size for monthly
rcParams['figure.figsize'] = 15, 5
plt.plot(date21, ter1, linestyle='-', linewidth=1.2, c='k', label=' MEGAN-Offline Monoterpene Emissions')
plt.plot(date21, v1, linestyle='--', linewidth=2.0, c='c', label=' MEGAN-Offline Isoprene Emissions')

# Create standard devation fill
plt.fill_between(date21, ter1-stdevjmeg2, ter1+stdevjmeg2, alpha=0.3, edgecolor='black', 
                 facecolor='black', linewidth=0.3 )

plt.fill_between(date21, v1-stdevjmeg1, v1+stdevjmeg1, alpha=0.3, edgecolor='black', 
                 facecolor='cyan', linewidth=0.3 )




#plt.plot(date21, rv1, linestyle='-', linewidth=3.0, c='b', label=' CTM Biogenic Emissions')
 #KOEH Data
#plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temperature')
#plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temperature')

# asthetics
plt.xlabel('Day')
plt.ylabel('Average Emissions (kg/$km^2$/hr)')
plt.title('Speciated Emissions February 2011')
plt.ylim(0, 25)
plt.xlim(0, 27)


# ticks
plt.xticks(range(1, 28, 1), [str(i) for i in range(1, 28, 1)])
plt.yticks(np.arange(0, 48, 2))  # use floats for y 

# Legend
plt.legend(loc='upper right')

# change plot size for daily
rcParams['figure.figsize'] = 7, 4

plt.show()
