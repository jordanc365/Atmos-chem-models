# -*- coding: utf-8 -*-
"""
Create emissions Timeseries MEGAN coupled to CTM

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

st = fm.variables['temp_a'][0:27, 14:24, 0, :, :]
st2 = fm.variables['temp_a'][0:27, 0:14, 0, :, :]
sg = fm.variables['skin_temp'][:, :, :, :]

# concatenate back to 24 hour period
vcon = np.concatenate((sv, sv2), axis=1)
tercon = np.concatenate((sm, sm2), axis=1)

tempcon = np.concatenate((st, st2), axis=1)

# Add emissions to find percentage isoprene and terpenes
totalf = (np.sum(vcon, axis=(2, 3)))
totalf1 = (np.sum(tercon, axis=(2, 3)))

# Mean over array
sv1 = vcon.mean(axis=(2, 3))
sm1 = tercon.mean(axis=(2, 3))
st1 = tempcon.mean(axis=(2, 3))
sg1 = sg.mean(axis=(2, 3))

# Reshape array for month
sr2 = np.reshape(sv1, 648)
smr2 = np.reshape(sm1, 648)
sr3 = np.reshape(st1, 648)
#sr4 = np.reshape(sg1, 675)

# Add monoterpenes and isoprene to same array Monthly
tot = [sum(x) for x in zip(sr2, smr2)]
totarr = np.array(tot)/9  # milimoles to kilograms

# Add monoterpenes and isoprene to same array daily
totdf = [sum(x) for x in zip(sv1, sm1)]
totarrdf = (np.array(totdf))/9  # milimoles to kilogram


# X axis dates instead of times
date = np.arange(totarr.shape[0])+1  # assume that delta time between data is 1
date31 = (date/24.)  # use days instead of hours

# change plot size
rcParams['figure.figsize'] = 15, 5

plt.plot(date31, sr3, linestyle='-', linewidth=3.0,
         c='r', label=' CTM inline MEGAN')

#plt.plot(ar, biosr2, linestyle='-', linewidth=3.0,
         #c='b', label=' CTM')  # KOEH Data

#plt.plot(date21, rv1, linestyle='-', linewidth=3.0,
         #c='c', label=' MEGAN')   # JMEGAN

# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temp')
# plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temperature')

# asthetics
plt.xlabel('Day')
plt.ylabel('Average Biogenic Emissions (kg/$km^2$/Hour)')
plt.title('Average Daily  Biogenic Emissions February 2011')
plt.ylim(14, 36)
plt.xlim(0, 27)


# ticks
plt.xticks(range(0, 28, 1), [str(i) for i in range(0, 28, 1)])
plt.yticks(np.arange(min(sr3)-1.75, max(sr3)-5, 1))  # use floats for y 

# Legend
#plt.legend(['CTM inline MEGAN', 'CTM', 'MEGAN '], loc='upper right', fontsize=10)

# change plot size
rcParams['figure.figsize'] = 7, 4

plt.show()
