# -*- coding: utf-8 -*-
"""
Create temp. Timeseries MEGAN coupled to CTM

@author: Jordan Capnerhurst 2016
"""

# import netCDF
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from pylab import rcParams
fm = netCDF4.Dataset('O:/Honours_data/KMEGAN/KFebmegv2.nc', 'r')

# plot Daily average
sv = (fm.variables['store_Megan'][0:27, 14:24, 0, 3:44, 1:43])  # Isoprene
sv2 = (fm.variables['store_Megan'][0:27, 0:14, 0, 3:44, 1:43])  # Isoprene
sm = (fm.variables['store_Megan'][0:27, 14:24, 1, 3:44, 1:43])  # Monoterp
sm2 = (fm.variables['store_Megan'][0:27, 0:14, 1, 3:44, 1:43])  # Monoterp

st = fm.variables['temp_a'][0:27, 14:24, 0, 3:44, 1:43]
st2 = fm.variables['temp_a'][0:27, 0:14, 0, 3:44, 1:43]
sg = fm.variables['skin_temp'][:, :, :, :]

# concatenate back to 24 hour period

tempcon = np.concatenate((st, st2), axis=1)

# reshape for month
tepresh = tempcon.reshape(648, 41, 42)

#calculate standard devation along lat and lons
stdevkmeg1 = np.std(tepresh, axis=(1, 2))

'''
# Add emissions to find percentage isoprene and terpenes
totalf = (np.sum(vcon, axis=(2, 3)))
totalf1 = (np.sum(tercon, axis=(2, 3)))
'''

# Mean over array

st1 = tepresh.mean(axis=(1, 2))



# X axis dates instead of times
date = np.arange(st1.shape[0])+1  # assume that delta time between data is 1
date31 = (date/24.)  # use days instead of hours

# change plot size
rcParams['figure.figsize'] = 15, 5
##############################################################################
plt.plot(date31, st1, linestyle='-', linewidth=2.0,
         c='r', label=' CSIRO-CTM-MEGAN')
         
# Create standard devation fill
plt.fill_between(date31, st1-stdevkmeg1, st1+stdevkmeg1, alpha=0.2, edgecolor='black', 
                 facecolor='red', linewidth=0.3)
##############################################################################

plt.plot(ar, t1, linestyle='-', linewidth=2.0,
         c='b', label=' CSIRO-CTM-Original')  # KOEH Data
         
# Create standard devation fill
plt.fill_between(ar, t1-stdevkoeh2, t1+stdevkoeh2, alpha=0.2, edgecolor='black', 
                 facecolor='blue', linewidth=0.3)

##############################################################################
plt.plot(date21, ter1, linestyle='-', linewidth=2.0,
         c='c', label=' MEGAN-Offline')   # JMEGAN
         
# Create standard devation fill
plt.fill_between(date21, ter1-stdev2, ter1+stdev2, alpha=0.25, edgecolor='black', 
                 facecolor='cyan', linewidth=0.3)
##############################################################################
                 
# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temp')
# plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temperature')

# asthetics
plt.xlabel('Day')
plt.ylabel('Temperature ($^o$c)')
plt.title('Average Daily  Biogenic Emissions February 2011')
plt.ylim(14, 44)
plt.xlim(0, 27)


# ticks
plt.xticks(range(0, 28, 1), [str(i) for i in range(0, 28, 1)])
plt.yticks(np.arange(min(st1)-1.06, max(st1)+2.94, 2))  # use floats for y 

# Legend
plt.legend( loc='upper right', fontsize=10)

# change plot size
rcParams['figure.figsize'] = 7, 4

plt.show()
