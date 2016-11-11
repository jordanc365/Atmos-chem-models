# -*- coding: utf-8 -*-
"""
Create Basemap without overlay

@author: Jordan Capnerhurst 2016
"""

# import things
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# Map
#  for 'low', not a numeral 1
map = Basemap(projection='merc', lat_0=-33, lon_0=151,
              resolution='h', area_thresh=0.1,
              llcrnrlon=147.804, llcrnrlat=-36.7246,  # Lower left corner
              urcrnrlon=153.114, urcrnrlat=-31.4146)  # Upper Right corner

plt.figure(figsize=(10, 10))

map.drawstates(color='black')
map.drawcoastlines()
map.drawcountries()
map.fillcontinents('white')
map.drawmapboundary()
map.drawrivers(color='lightblue')
map.shadedrelief()

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 3), labels=[0, 0, 0, 1], fontsize=10)
map.drawparallels(np.arange(-90, 90, 3), labels=[1, 1, 1, 1], fontsize=10)

# City markers and names
lons = [151.2070, 150.8931, 151.7789, 149.1287]
lats = [-33.8675, -34.4250, -32.9267, -35.2820]
x, y = map(lons, lats)
map.plot(x, y, 'ko', markersize=15)

labels = ['      Sydney', 'Wollongong', '  Newcastle', '   Canberra']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-90000, ypt+10000, label)
fontsize = 14

# North Arrow
lons = [152.7]
lats = [-36.5]
x, y = map(lons, lats)
map.plot(x, y, 'k^', markersize=30)

# Add scale bar
map.drawmapscale(151, -36.4, 151, -34, 100, barstyle='simple', units='km',
                 fontsize=15, yoffset=None, labelstyle='simple', fontcolor='k',
                 fillcolor1='k', fillcolor2='k', ax=None, format='%d',
                 zorder=None)

labels = ['North', '', '  ', '   ']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-20000, ypt+30000, label)
fontsize = 40

plt.title("Sydney Metropolitan Region", fontsize=20)
