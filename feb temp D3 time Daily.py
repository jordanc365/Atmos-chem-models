# -*- coding: utf-8 -*-
"""
Create temperature Timeseries CTM 2013

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

totalz = (np.sum(vcon, axis=(2, 3)))
# Mean over array
bios1 = totalz.mean(axis=(0))
t1 = tempcon.mean(axis=(2, 3))
g1 = g.mean(axis=(2, 3))

stdevkoeh = np.std(totalz, axis=(0))

# Reshape array for month
#biosr2 = np.reshape(t1, 648)  # Due to Dimensions being 3x3
#r3 = np.reshape(t1, 648)
#r4 = np.reshape(g1, 675)


# X axis dates instead of times for monthly 
#date = np.arange(biosr2.shape[0])  # assume  delta time between data is 1 hour
#date1 = (date/24.)  # use days instead of hours
#ar = np.array(date1)
time= np.arange(0, 24, 1)

# change plot size for monthly
#rcParams['figure.figsize'] = 15, 5
  

#############################################################################################                    
plt.plot(time, no, linestyle='-', linewidth=2.0, c='b',
         label='February 2011')

# Create standard devation fill
plt.fill_between(time,  no-(stdevkoeh2*0.5),  no+(stdevkoeh2*0.5), alpha=0.2, edgecolor='black', 
                 facecolor='blue', linewidth=0.3)
#############################################################################################            
plt.plot(time, bios2, linestyle='-', linewidth=2.0, c='r',
         label='February 2013')

# Create standard devation fill
plt.fill_between(time,  bios2-stdevkoeh2,  bios2+stdevkoeh2, alpha=0.2, edgecolor='black', 
                 facecolor='red', linewidth=0.3)
                 
#############################################################################################                  
# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temp')
# plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temp')

# asthetics
plt.xlabel('Hour (UTC +10')
plt.ylabel('Average Biogenic Emissions (kg/hour)')
#plt.title('Monthly Ambiant temperature Feburary 2011 from CTM')
plt.ylim(0, 150000)
plt.xlim(0, 23)


# tick
plt.xticks(range(0, 24, 1), [str(i) for i in range(0, 24, 1)])
plt.yticks(np.arange(0, 210000,10000))  # use floats for y value

# Legend
plt.legend(loc='upper left')

# change plot size
rcParams['figure.figsize'] = 7, 4

plt.show()
# plt.savefig('./bmap_syd.png')
