# -*- coding: utf-8 -*-
"""
Create average emissions map 2013

@author: Jordan Capnerhurst 2016
"""

# import things
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from mpl_toolkits.basemap import maskoceans

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
# [ time, source, lon, lat ] Feb
# plot Daily average
m21 = m2.variables['store_Bio'][:, 14:24, 0, 22:36, 25:45]
m22 = m2.variables['store_Bio'][:, 0:14, 0, 22:36, 25:45]

# concatenate back to 24 hour period
feb = np.concatenate((m21, m22), axis=1)/81

# average
febav = np.mean(feb, axis=(0, 1))

# Reshape for month
#febres = feb.reshape(672, 60, 60) 



####################################################
# [ time, source, lon, lat ] Feb
# plot Daily average
mt21 = m2.variables['temp_a'][:, 14:24, 0, 22:36, 25:45]
mt22 = m2.variables['temp_a'][:, 0:14, 0, 22:36, 25:45]

# concatenate back to 24 hour period
febt = np.concatenate((mt21, mt22), axis=1)

# average
febtav = np.mean(febt, axis=(0, 1))

# Reshape for month
#febres = feb.reshape(672, 60, 60) 

datalats = m3.variables['lat'][:] 
datalons = m2.variables['lon'][:]

cc = m1.variables['cloud_cover'][0:28, 0:23, 0, :, :]
cc1 = cc.mean(axis=(0, 1))*100


##############################################################################
# Define map parameters

bmax = 9
bmin = 0
trans = 1

##############################################################################

# Other variables
at = (m2.variables['lndtype'][0, :, :])
st = (m2.variables['soiltype'][0, :, :])
lai = (m1.variables['lai'][0, :, :])
t = (m2.variables['temp_a'][:, 0:24, 0, :, :])
g = (m2.variables['skin_temp'][:, 0:24, :, :])


'''
# Mean over array
v1 = m2con.mean(axis=(0, 1))
t1 = t.mean(axis=(0, 1))
g1 = g.mean(axis=(0, 1))
#cc1 = cc.mean(axis=(0, 1))
#cr1 = cr.mean(axis=(0, 1))
#ft1 = cr.mean(axis=(0, 1))
'''



# Map
#  for 'low', not a numeral 1
map = Basemap(projection='merc', lat_0=-33, lon_0=151,
              resolution='h', area_thresh=0.1,
              llcrnrlon=147.804, llcrnrlat=-36.7246,  # Lower left corner
              urcrnrlon=153.114, urcrnrlat=-31.4146)  # Upper Right corner


plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# Mask oceans and dams
#mocedata = maskoceans(mlons, mlats, febav,  inlands=True, resolution='h',
                      #grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats,cc1, latlon=True, zorder=1, vmin=bmin, vmax=bmax,
               cmap=cool, alpha=trans)


map.drawstates(color='white', linewidth=3)
map.drawcoastlines(color='white', linewidth=3)
# map.drawcountries()
# map.fillcontinents('white')
map.drawmapboundary()
# map.drawrivers(color='black', linewidth=2)
# map.shadedrelief()


# City markers and names
lons = [151.2070, 150.8931, 151.7789, 149.1287]
lats = [-33.8675, -34.4250, -32.9267, -35.2820]
x, y = map(lons, lats)
map.plot(x, y, 'wo', markersize=15)

labels = ['   Sydney', '   Wollongong', ' Newcastle', '    Canberra']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt+5000, ypt-9000, label, color='white', fontsize=20)

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 1), labels=[0, 0, 0, 1], fontsize=10,
                  color='white', linewidth=2)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=10,
                  color='white', linewidth=2)

# North Arrow
lons = [152.7]
lats = [-36.5]
x, y = map(lons, lats)
map.plot(x, y, 'w^', markersize=30)

# Add scale bar
map.drawmapscale(151, -36.4, 151, -34, 50, barstyle='simple', units='km',
                 fontsize=15, yoffset=None, labelstyle='simple', fontcolor='w',
                 fillcolor1='w', fillcolor2='w', ax=None, format='%d',
                 zorder=None)

labels = ['North', '', '  ', '   ']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-20000, ypt+30000, label, color='white')
fontsize = 40


# Add colour bar
col = map.pcolormesh(mlons, mlats, cc1,
                     latlon=True, zorder=1, vmin=bmin, vmax=bmax, 
                     cmap=cool, alpha=trans)  

cb = map.colorbar(col, "bottom", size="5%", pad="2%")
cb.set_label('Monthly Average Ambient Temperature ($^o$c)', fontsize=17)

plt.title("CTM Domain 3 3x3 Leaf Area Index\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)

#map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)

plt.show()
# plt.savefig('./bmap_syd.png')
