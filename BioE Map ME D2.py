# -*- coding: utf-8 -*-
"""
Create Map of emissions MEGAN 3x3 grid

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

f = netCDF4.Dataset('O:/Honours_data/JMEGAN/FebJMEGv2.nc', 'r')

dim = netCDF4.Dataset('O:/Honours_data/JMEGAN/GRIDCRO2D_d03.nc',
                      'r')
                      
# [ time, source, lon, lat ]
v = (f.variables['ISOP'][1:29, 15:24, :, :]*68.12)
v2 = (f.variables['ISOP'][1:29, 0:15, :, :]*68.12)
ter = (f.variables['TERP'][1:29, 15:24, :, :]*136.298)
ter2 = (f.variables['TERP'][1:29, 0:15, :, :]*136.298)

# concatenate back to 24 hour period
vcon = np.concatenate((v, v2), axis=1)
tercon = np.concatenate((ter, ter2), axis=1)
#Etho = (f.variables['ETOH'][:, :, 0, :, :])*46.068 # Ethonol
#Metho = (f.variables['MEOH'][:, :, 0, :, :])*32.04 # Methanol
#Alde = (f.variables['ALDX'][:, :, 0, :, :])*44.052 # Aldehyde
#Ace = (f.variables['ALD2'][:, :, 0, :, :])*58.08 # Acetone

vhalf= vcon[:, 6:20, :, :]
#vhalf1= vcon[:, 20:23, :, :]

vhalf2= tercon[:, 6:20, :, :]
#vhalf21= tercon[:, 20:23, :, :]

#vcon2 = np.concatenate((vhalf, vhalf1), axis=1)
#tercon2 = np.concatenate((vhalf2, vhalf21), axis=1)




# Mean over array
v1 = vhalf.mean(axis=(0, 1))
ter1 = vhalf2.mean(axis=(0, 1))
#t1 = t.mean(axis=(0, 1))
#g1 = g.mean(axis=(0, 1))
#cc1 = cc.mean(axis=(0, 1))
#cr1 = cr.mean(axis=(0, 1))
#ft1 = cr.mean(axis=(0, 1))

# Pearson R stats via for loop
# r1_array = np.zeros((60, 60))
# for v1 in range(0, 61):
  #  print(v1)
   # for t1 in range(0, 61):
       # f = 'cheese'

# Add emissions to get total emissions
#total = np.sum(v, axis=(0, 1, 2, 3))*3.6
#total1 = np.sum(ter, axis=(0, 1, 2, 3))*3.6

# Add monoterpenes and isoprene to same array daily
totdf = [sum(x) for x in zip(v1, ter1)]
totarrdf = np.array(totdf)*3.6/9 # milimoles to kilograms

datalats = f.variables['lat'][:]
datalons = f.variables['lon'][:]

# Map
#  for 'low', not a numeral 1
map = Basemap(projection='merc', lat_0=-33, lon_0=151,
              resolution='h', area_thresh=0.1,
              llcrnrlon=150.079, llcrnrlat=-34.836,  # Lower left corner
              urcrnrlon=151.849, urcrnrlat=-33.066)  # Upper Right corner


plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, totarrdf,  inlands=True, resolution='h',
                      grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, latlon=True, vmax=15,
                zorder=1, cmap=cool)
# alpha=0.8)

map.drawstates(color='black', linewidth=3)
map.drawcoastlines(color='black', linewidth=3)

# map.drawcountries()
# map.fillcontinents('white')
map.drawmapboundary()
# map.drawrivers(color='black', linewidth=2)
# map.shadedrelief()

# City markers and names
lons = [151.2070, 150.8931, 151.7789, 149.1287]
lats = [-33.8675, -34.4250, -32.9267, -35.2820]
x, y = map(lons, lats)
map.plot(x, y, 'ro', markersize=15)

labels = ['   Sydney', '   Wollongong', '   ', '']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt+5000, ypt-9000, label, color='black', fontsize=20)

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 1), labels=[0, 0, 0, 1], fontsize=10,
                  color='black', linewidth=2)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=10,
                  color='black', linewidth=2)

# North Arrow
lons = [151.3]
lats = [-34.67]
x, y = map(lons, lats)
map.plot(x, y, 'k^', markersize=30)

labels = ['North', '', '  ', '   ']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-10000, ypt+8000, label, color='black', fontsize=20)

# Add scale bar
map.drawmapscale(151.5, -34.65, 151, -34, 5, barstyle='simple', units='km',
                 fontsize=15, yoffset=None, labelstyle='simple', fontcolor='k',
                 fillcolor1='k', fillcolor2='k', ax=None, format='%d',
                 zorder=None)

# Add colour bar
col = map.pcolormesh(mlons, mlats, mocedata, latlon=True,cmap=cool, vmax=15)

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Average Total Biogenic Emissions (kg/$km^2$/hr)', fontsize=15)

plt.title("Megan Domain 3 3x3  Monthly Distribution of Total Biogenic \
Emissions \n Sydney Metropolitan Region Feburary 2011", fontsize=10)


# plt.savefig('bmap_syd.png', dpi=100, bbox_inches='tight')
plt.show()

