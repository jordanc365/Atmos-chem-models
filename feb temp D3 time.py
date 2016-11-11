# -*- coding: utf-8 -*-
"""
Create D3 monthly emissions timeseries CTM 2013

@author: Jordan Capnerhurst 2016
"""

# import things
import matplotlib.pyplot as plt
import numpy as np
import netCDF4


# Set colour map
cool = cm = plt.get_cmap('jet')

# Import data
m1 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JAN.nc', 'r')  # January
m2 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/FEB.nc', 'r')  # Febuary
m3 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/MAR.nc', 'r')  # March
m4 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/APR.nc', 'r')  # April
m5 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/MAY.nc', 'r')  # May
m6 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JUN.nc', 'r')  # June
m7 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JUL.nc', 'r')  # July
m8 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/AUG.nc', 'r')  # August
m9 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/SEP.nc', 'r')  # September
m10 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/OCT.nc', 'r')  # October
mnov = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/NOV.nc', 'r')  # November
mdec = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/DEC.nc', 'r')  # Decmber


# Import variables
# Store_bio

##############################################################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011
# plot Daily average
m21 = m2.variables['store_Bio'][:27, 14:24, 0, 22:36, 25:45]
m22 = m2.variables['store_Bio'][:27, 0:14, 0, 22:36, 25:45]

# concatenate back to 24 hour period
feb = np.concatenate((m21, m22), axis=1)

# average
febav = np.mean(feb, axis=(0, 1))

#febavres = np.reshape(febav, 3600)

# Reshape for month
febres = feb.reshape(648, 14, 20) 

stdfeb = np.std(febres, axis=(1, 2))

febresav = np.mean(febres, axis=(1,2))

# Daily average

totalz2 = (np.sum(feb, axis=(2, 3)))
bios2 = totalz2.mean(axis=(0))
stdevkoeh2 = np.std(totalz2, axis=(0))*2

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt21 = m2.variables['temp_a'][:, 14:24, 0, 22:36, 25:45]
mt22 = m2.variables['temp_a'][:, 0:14, 0, 22:36, 25:45]

# concatenate back to 24 hour period
febt = np.concatenate((mt21, mt22), axis=1)

# average
febtav = np.mean(febt, axis=(0, 1))

#febtavres = np.reshape(febtav, 3600)

# Reshape for month
febtres = febt.reshape(672, 14, 20) 

stdfebt = np.std(febtres, axis=(1, 2))

febtresav = np.mean(febtres, axis=(1,2))


##############################################################################
# Define map parameters

max1 = 9
min1 = 0
trans = 0.5

##############################################################################

# Other variables
at = (m2.variables['lndtype'][0, :, :])
st = (m2.variables['soiltype'][0, :, :])
lai = (m1.variables['lai'][0, :, :])
laires = at.reshape(3600)
t = (m2.variables['temp_a'][:, 0:24, 0, :, :])
g = (m2.variables['skin_temp'][:, 0:24, :, :])

##############################################################################

# X axis dates instead of times for monthly 
a = np.arange(febtresav.shape[0])  # assume  delta time between data is 1 hour
date1 = (date/24.)  # use days instead of hours
ar = np.array(date1)

plt.figure(figsize=(15, 5))

##############################################################################

plt.plot(ar, febtresav, linestyle='-', linewidth=2, c='r',
         label='Average Ambient Temperature ($^o$c)')

# Create standard devation fill
plt.fill_between(ar,  febtresav-stdfebt,  febtresav+stdfebt, alpha=0.3, edgecolor='black', 
                 facecolor='red', linewidth=0.3)
        
############################################################################################################################################################                                 

plt.plot(ar, febresav, linestyle='-', linewidth=2, c='b',
         label='Average Total Biogenic Emissions (kg/$km^2$/hr)')

# Create standard devation fill
plt.fill_between(ar,  febresav-stdfeb,  febresav+stdfeb, alpha=0.2, edgecolor='black', 
                 facecolor='blue', linewidth=0.3)

############################################################################################################################################################
plt.title("CTM Domain 3 3x3 Leaf Area Index\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)

plt.xlabel('Month')
#plt.ylabel('Average Ambient Temperature ($^o$c)')

plt.ylim(0, 48)
plt.xlim(0, 27)

# tick
plt.xticks(range(0, 28, 1 ), [str(i) for i in range(0, 28, 1)])
plt.yticks(np.arange(0, 49, 2))  # use floats for y value

#plt.legend(loc='upper center', fontsize=10)

plt.show()
# plt.savefig('./bmap_syd.png')
