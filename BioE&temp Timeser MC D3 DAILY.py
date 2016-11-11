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

# standardss devitations
stdevkmeg1 = np.std(totalf, axis=(0))
stdevkmeg2 = np.std(totalf1, axis=(0))

# Mean over array
sv1 = totalf.mean(axis=(0))
sm1 = totalf1.mean(axis=(0))
st1 = tempcon.mean(axis=(2, 3))
sg1 = sg.mean(axis=(2, 3))

hour = np.arange(0, 24, 1)

# Add monoterpenes and isoprene to same array daily
totdf = [sum(x) for x in zip(sv1, sm1)]
totarrdf = (np.array(totdf))  # milimoles to kilogram

# stds for std devs
tot2 = [sum(x) for x in zip(stdevkmeg1, stdevkmeg2)]
totstdkmeg = np.array(tot2)  

##############################################################################
'''
plt.plot(totarrdf, linestyle='-', linewidth=2.0,
         c='r', label=' CSIRO-CTM-MEGAN')
         
# Create standard devation fill
plt.fill_between(hour, totarrdf-totstdkmeg, totarrdf+totstdkmeg, alpha=0.5, edgecolor='black', 
                 facecolor='red', linewidth=0.3)
'''
##############################################################################
plt.plot(bios1, linestyle='-', linewidth=2.0,
         c='b', label=' February 2011')  # KOEH Data
         
# Create standard devation fill
plt.fill_between(hour, bios1-stdevkoeh, bios1+stdevkoeh, alpha=0.3, edgecolor='black', 
                 facecolor='blue', linewidth=0.3)

##############################################################################
'''
plt.plot(totarrd, linestyle='-', linewidth=2.0,
         c='c', label=' MEGAN-Offline')   # JMEGAN
         
# Create standard devation fill
plt.fill_between(hour, totarrd-stdtot, totarrdf+stdtot, alpha=0.3, edgecolor='black', 
                 facecolor='cyan', linewidth=0.3)
'''              
############################################################################################################################################################
plt.plot(bios2, linestyle='-', linewidth=2.0,
         c='r', label=' February 2013')   # CTM Year 2013
         
# Create standard devation fill
plt.fill_between(hour, bios2-stdevkoeh2, bios2+stdevkoeh2, alpha=0.3, edgecolor='black', 
                 facecolor='red', linewidth=0.3)
            
############################################################################################################################################################

# plt.plot(t1, linestyle='--', linewidth=5.0, c='r', label='Ambiant Temp')
# plt.plot(g1, linestyle='--', linewidth=5.0, c='b', label='Skin Temperature')

# asthetics
plt.xlabel('Hour (UTC +10)')
plt.ylabel('Average Biogenic Emissions (kg/Hour)')
#plt.title('Average Daily  Biogenic Emissions February 2011')
plt.ylim(0, 220000)
plt.xlim(0, 23)


# ticks
plt.xticks(range(0, 24, 1), [str(i) for i in range(0, 24, 1)])
plt.yticks(np.arange(min(totarrdf)-1352, max(totarrdf)+120000, 10000))  # use floats for y 

# Legend
plt.legend( loc='upper left', fontsize=8)

# change plot size
rcParams['figure.figsize'] = 7, 4

plt.show()
