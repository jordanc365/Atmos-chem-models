# -*- coding: utf-8 -*-
"""
Create temp. and emissions Timeseries MEGAN coupled to CTM

@author: Jordan Capnerhurst 2016
"""

# import netCDF
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from pylab import rcParams
fm = netCDF4.Dataset('O:/Honours_data/KMEGAN/KFebmegv2.nc', 'r')

# plot Daily average
sv = (fm.variables['store_Megan'][1:27, 14:24, 0, 3:44, 1:43])  # Isoprene
sv2 = (fm.variables['store_Megan'][1:27, 0:14, 0, 3:44, 1:43])  # Isoprene
sm = (fm.variables['store_Megan'][1:27, 14:24, 1, 3:44, 1:43])  # Monoterp
sm2 = (fm.variables['store_Megan'][1:27, 0:14, 1, 3:44, 1:43])  # Monoterp

st = fm.variables['temp_a'][1:27, 14:24, 0, 3:44, 1:43]
st2 = fm.variables['temp_a'][1:27, 0:14, 0, 3:44, 1:43]

# concatenate back to 24 hour period
vcon = np.concatenate((sv, sv2), axis=1)
tercon = np.concatenate((sm, sm2), axis=1)
tempcon = np.concatenate((st, st2), axis=1)

# Reshape array for month
sr2 = vcon.reshape(624, 41, 42)
smr2 = tercon.reshape(624, 41, 42)
sr3 = tempcon.reshape(624, 41, 42)

#calculate standard devation along lat and lons
stdevkmeg1 = np.std(sr2, axis=(1, 2))
stdevkmeg2 = np.std(smr2, axis=(1, 2))
stdevkmeg3 = np.std(sr3, axis=(1, 2))


# Add emissions to find percentage isoprene and terpenes
#totalf = (np.sum(vcon, axis=(2, 3)))
#totalf1 = (np.sum(tercon, axis=(2, 3)))

# Mean over array
sv1 = sr2.mean(axis=(1, 2))/9
sm1 = smr2.mean(axis=(1, 2))/9
st1 = sr3.mean(axis=(1, 2))


# Add monoterpenes and isoprene to same array Monthly
tot = [sum(x) for x in zip(sv1, sm1)]
totarr = np.array(tot)  # milimoles to kilograms

# Add monoterpenes and isoprene to same array Monthly
tot = [sum(x) for x in zip(sv1, sm1)]
totarr = np.array(tot)/9  # milimoles to kilograms

# stds for total emissions
totdf = [sum(x) for x in zip(stdevkmeg1, stdevkmeg2)]
totstdkmeg = (np.array(totdf))/9 # milimoles to kilogram


# X axis dates instead of times
date = np.arange(sm1.shape[0])+1  # assume that delta time between data is 1
date31 = (date/24.)  # use days instead of hours

# change plot size
rcParams['figure.figsize'] = 15, 5

##############################################################################
plt.plot(date31, sm1, linestyle='-', linewidth=2.0,
         c='r', label=' CSIRO-CTM-MEGAN')
         
# Create standard devation fill
plt.fill_between(date31, totarr-totstdkmeg, totarr+totstdkmeg, alpha=0.3, edgecolor='black', 
                 facecolor='red', linewidth=0.3)
                 
##############################################################################
plt.plot(date31, bios1, linestyle='-', linewidth=2.0,
         c='b', label='CSIRO-CTM-Original')  # KOEH Data
         
# Create standard devation fill
#plt.fill_between(ar, bios1-stdevkoeh1, bios1+stdevkoeh1, alpha=0.3, edgecolor='black', 
                # facecolor='blue', linewidth=0.3)
         
##############################################################################
plt.plot(date21, totarr2, linestyle='-', linewidth=2.0,
         c='c', label='MEGAN-Offline')   # JMEGAN
         
# Create standard devation fill
plt.fill_between(date21, totarr2-totstdjmeg, totarr2+totstdjmeg, alpha=0.3, edgecolor='black', 
                 facecolor='cyan', linewidth=0.3)
                 
##############################################################################
# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temp')
# plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temperature')

# asthetics
plt.xlabel('Day')
plt.ylabel('Average Total Biogenic Emissions (kg/$km^2$/hr)')
plt.title('Average Daily  Biogenic Emissions February 2011')
plt.ylim(0, 48)
plt.xlim(0, 27)


# ticks
plt.xticks(range(0, 28, 1), [str(i) for i in range(0, 28, 1)])
plt.yticks(np.arange(min(totarr2)-0.04, max(totarr2)+25.94, 2))  # use floats for y 

# Legend
plt.legend(loc='upper right', fontsize=10)

# change plot size
rcParams['figure.figsize'] = 25, 5

plt.show()
