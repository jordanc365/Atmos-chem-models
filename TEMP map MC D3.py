# -*- coding: utf-8 -*-
"""
Create temperature map MEGAN coupled to CTM

@author: Jordan Capnerhurst 2016
"""

# import things
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4


# Set colour map
cool = cm = plt.get_cmap('jet')

# Leaf area index used in model
#flai = netCDF4.Dataset('O:/Honours_data/KMEGAN/syd_3km_laifeb.nc', 'r')

f = netCDF4.Dataset('O:/Honours_data\KMEGAN\FebMEG.nc',
                    'r')

#at = (f.variables['lndtype'][0, :, :])
#st = (f.variables['soiltype'][0, :, :])
t = f.variables['temp_a'][1:, 0:23, 0, 3:44, 1:55]*1.025
g = f.variables['skin_temp'][1:, 0:23, 3:44, 1:55]


# Mean over array
t1 = t.mean(axis=(0, 1))
g1 = g.mean(axis=(0, 1))

datalats = f.variables['lat'][3:44]
datalons = f.variables['lon'][1:55]

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
map.pcolormesh(mlons, mlats, t1,  latlon=True, zorder=1, vmax=24.8, vmin=18.4,
               cmap=cool, alpha=trans)

'''
map.drawstates(color='black', linewidth=3)
map.drawcoastlines(color='black', linewidth=3)
# map.drawcountries()
# map.fillcontinents('white')
map.drawmapboundary()
map.drawrivers(color='red', linewidth=10, zorder=10)
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
'''
# Add colour bar
col = map.pcolormesh(mlons, mlats, t1, latlon=True, vmax=24.8, vmin=18.4,
                     cmap=cool, alpha=trans)
'''
# Add scale bar
map.drawmapscale(151.5, -34.65, 151, -34, 5, barstyle='simple', units='km',
                 fontsize=15, yoffset=None, labelstyle='simple', fontcolor='k',
                 fillcolor1='k', fillcolor2='k', ax=None, format='%d',
                 zorder=None)
'''

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Average Total Biogenic Emissions (kg/$km^2$/hr)', fontsize=15)

plt.title("CTM Inline MEGAN Domain 3 3x3 Monthly Distribution of \
Biogenic Emissions \n Sydney Metropolitan Region  February 2011", fontsize=10)

map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)
plt.show()
# plt.savefig('./bmap_syd.png')


