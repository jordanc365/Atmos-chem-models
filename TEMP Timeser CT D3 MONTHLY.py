# -*- coding: utf-8 -*-
"""

Create Temperature Timeseries CTM

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
t = fc.variables['temp_a'][:, 14:24, 0, 6:46, 6:48]
t2 = fc.variables['temp_a'][:, 0:14, 0, 6:46, 6:48]


# concatenate back to 24 hour period
vcon = np.concatenate((vbios, vbios2), axis=1)
tempcon = np.concatenate((t, t2), axis=1)

# Reshape array for month
biosr1 = vcon.reshape(648, 40, 42)  
r3 = tempcon.reshape(648, 40, 42)

#calculate standard devation along lat and lons
stdevkoeh1 = np.std(biosr1, axis=(1, 2))/9
stdevkoeh2 = np.std(r3, axis=(1, 2))

# add emissions over arrays
#totalz = (np.sum(vcon, axis=(2, 3)))*0.91

# Mean over array
bios1 = biosr1.mean(axis=(1, 2))/9
t1 = r3.mean(axis=(1, 2))



# X axis dates instead of times for monthly 
date = np.arange(bios1.shape[0])  # assume  delta time between data is 1 hour
date1 = (date/24.)  # use days instead of hours
ar = np.array(date1)

# change plot size for monthly
rcParams['figure.figsize'] = 25, 5

plt.plot(ar, t1, linestyle='-', linewidth=2.0, c='b',
         label='Biogenic Emissions')
         
# Create standard devation fill
plt.fill_between(ar, t1-stdevkoeh2, t1+stdevkoeh2, alpha=0.3, edgecolor='black', 
                 facecolor='blue', linewidth=0.3, label='1 Standard devation' )
         

# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temp')
# plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temp')

# asthetics
plt.xlabel('Day')
plt.ylabel('Temperature ($^o$c)')
plt.title('Monthly Ambiant temperature Feburary 2011 from CTM')
#plt.ylim(14, 26)
plt.xlim(0, 27)


# tick
plt.xticks(range(0, 28, 1), [str(i) for i in range(0, 28, 1)])
#plt.yticks(np.arange(min(r3)-1.45, max(r3)+2, 1))  # use floats for y value

# Legend
plt.legend(['Store Bio', 'Ambiant temp', 'Skin Temp'], loc='upper left')

# change plot size
rcParams['figure.figsize'] = 7, 4

plt.show()
# plt.savefig('./bmap_syd.png')
