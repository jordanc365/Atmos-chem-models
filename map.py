# -*- coding: utf-8 -*-
"""
Create Basic Basemap with pcolourmesh

@author: Jordan Capnerhurst 2016
"""

# import things
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4


f = netCDF4.Dataset('C:/Honours_data/gmr_20130201_cb05_aer2_0.090dg_feb13.nc',
                    'r')
# [ time, source, lon, lat ]
v = (f.variables['store_Bio'][0:24, 0, 0:60, 0:60])/81

datalons = f.variables['lon'][:]
datalats = f.variables['lat'][:]

# Map
#  for 'low', not a numeral 1
map = Basemap(projection='merc', lat_0=-33, lon_0=151,
              resolution='h', area_thresh=0.1,
              llcrnrlon=147.804, llcrnrlat=-36.7246,  # Lower left corner
              urcrnrlon=153.114, urcrnrlat=-31.4146)  # Upper Right corner

plt.figure(figsize=(10, 10))


mlons, mlats = np.meshgrid(datalons, datalats)

# Mean over array
v1 = v.mean(axis=(0))


map.pcolormesh(mlons, mlats, v1, vmin=0, vmax=5.6, latlon=True,
               zorder=1)  # , alpha=0.8)

map.drawstates(color='black', linewidth=3)
map.drawcoastlines(color='black', linewidth=3)
# map.drawcountries()
# map.fillcontinents('white')
map.drawmapboundary()
map.drawrivers(color='black', linewidth=2)
# map.shadedrelief()

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 3), labels=[0, 0, 0, 1], fontsize=10,
                  color='white', linewidth=5)
map.drawparallels(np.arange(-90, 90, 3), labels=[1, 0, 0, 0], fontsize=10,
                  color='white', linewidth=5)

# City markers and names
lons = [151.2070, 150.8931, 151.7789, 149.1287]
lats = [-33.8675, -34.4250, -32.9267, -35.2820]
x, y = map(lons, lats)
map.plot(x, y, 'wo', markersize=15)

labels = ['Sydney', 'Wollongong', 'Newcastle', 'Canberra']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt+100, ypt-35000, label, color='white', fontsize=20)


# North Arrow
lons = [152.6]
lats = [-36.5]
x, y = map(lons, lats)
map.plot(x, y, 'w^', markersize=30)

labels = ['North', '', '  ', '   ']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-40000, ypt+30000, label, color='white', fontsize=20)

# Add colour bar
col = map.pcolormesh(mlons, mlats, v1, latlon=True, zorder=1)
cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Biogenic emissions kg/hr/$kg^2$', fontsize=15)

plt.title("AverageBiogenic emissions \n Sydney Metropolitan Region 1/2/2013",
          fontsize=20)

plt.show()
# plt.savefig('./bmap_syd.png')
