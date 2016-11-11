# -*- coding: utf-8 -*-
"""
Create LAI map for CTM coupled to MEGAN

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

f = netCDF4.Dataset('O:/Honours_data/KMEGAN/syd_3km_laifeb.nc','r')

v1 = f.variables['lai'][:, :]
lait = np.array(v1)
datalons = f.variables['lon'][:]
datalats = f.variables['lat'][:]

# Map
#  for 'low', not a numeral 1
map = Basemap(llcrnrlon=150.079, llcrnrlat=-34.836,  # Lower left corner
              urcrnrlon=151.849, urcrnrlat=-33.066,  epsg=4269)

plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

trans = 0.2

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, lait,  inlands=True, resolution='h',
                      grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, latlon=True, zorder=1, cmap=cool,  vmin=0, vmax=6, alpha=trans)

'''
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
col = map.pcolormesh(mlons, mlats, mocedata, latlon=True, zorder=1, cmap=cool, vmin=0, vmax=6, alpha=trans)

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('LAI ($m^2$/$m^2$)', fontsize=15)

plt.title("MEGAN CTM Domain 3  3x3 Leaf Area Index \n Sydney Metropolitan Region Feburary 2011",
          fontsize=10)

map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)

plt.show()
# plt.savefig('./bmap_syd.png')
