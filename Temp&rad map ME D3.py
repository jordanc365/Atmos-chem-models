# -*- coding: utf-8 -*-
"""
Radation and temperature map MEGAN

@author: Jordan Capnerhurst 2016
"""

# import things
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4

# Set colour map
cool = cm = plt.get_cmap('jet')

f = netCDF4.Dataset('O:/Honours_data/JMEGAN/rad1.nc','r')

dim = netCDF4.Dataset('O:/Honours_data/JMEGAN\Grid\GRIDCRO2D_d03.nc',
                      'r')

# [ time, source, lon, lat ]
rad = (f.variables['GSW'][:, 14:24, 0, :, 10:155]) # Radiation at ground level
rad2 = (f.variables['GSW'][:, 0:14, 0, :, 10:155]) 
temp = (f.variables['TEMP2'][:, 14:24, 0, :, 10:155]) # temp at 2m
temp2 = (f.variables['TEMP2'][:, 0:14, 0, :, 10:155])

# concatenate back to 24 hour period
vcon = np.concatenate((rad, rad2), axis=1)
tercon = np.concatenate((temp, temp2), axis=1)


# Averages
v1 = vcon.mean(axis=(0, 1))
ter1 = tercon.mean(axis=(0 ,1))-273.15

datalats = f.variables['LAT'][0, 0, :, 0]
datalons = f.variables['LON'][0, 0, 0, 10:155]

# Map
#  for 'low', not a numeral 1
map = Basemap(llcrnrlon=150.1,
              llcrnrlat=-34.7177,
              urcrnrlon=151.651,
              urcrnrlat=-33.5651, epsg=4269)



plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

trans= 0.2

# Colour mesh
map.pcolormesh(mlons, mlats, ter1, latlon=True, vmax=24.8,
                zorder=1, cmap=cool, alpha=trans)
'''
map.drawstates(color='black', linewidth=3)
map.drawcoastlines(color='black', linewidth=3)

map.drawcountries()
# map.fillcontinents('white')
#map.drawmapboundary()
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
'''

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 1), labels=[0, 0, 0, 1], fontsize=10,
                  color='black', linewidth=2)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=10,
                  color='black', linewidth=2)
'''
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
'''

# Add colour bar
col = map.pcolormesh(mlons, mlats, ter1, latlon=True, cmap=cool, vmax=24.8, vmin=18.4, alpha=trans)

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Normalised Biogenic Emissions', fontsize=15)

plt.title("Megan Domain 3 1x1  Monthly  Total Biogenic \
Emissions \n Sydney Metropolitan Region Feburary 2011", fontsize=10)

map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)
# plt.savefig('bmap_syd.png', dpi=100, bbox_inches='tight')
plt.show()
