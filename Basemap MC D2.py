# -*- coding: utf-8 -*-
"""
Create Basemap without overlay MEGAN coupled to CTM Domain 2

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
              llcrnrlon=150.079, llcrnrlat=-34.836,  # Lower left corner
              urcrnrlon=151.849, urcrnrlat=-33.066)  # Upper Right corner

plt.figure(figsize=(10, 10))

map.drawstates(color='black')
map.drawcoastlines()
map.drawcountries()
map.fillcontinents('white')
map.drawmapboundary()
map.drawrivers(color='lightblue')
map.shadedrelief()

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, ), labels=[0, 0, 0, 1], fontsize=10)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 1, 1, 1], fontsize=10)

# City markers and names
lons = [151.2070, 150.8931, 151.7789, 149.1287]
lats = [-33.8675, -34.4250, -32.9267, -35.2820]
x, y = map(lons, lats)
map.plot(x, y, 'ko', markersize=15)

labels = ['      Sydney', 'Wollongong', '  ', '  New ']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt+10000, ypt+1000, label)
fontsize = 14

# North Arrow
lons = [151.72]
lats = [-34.7]
x, y = map(lons, lats)
map.plot(x, y, 'k^', markersize=30)

labels = ['North', '', '  ', '   ']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-11300, ypt-15000, label, color='black', fontsize=20)

# Add scale bar
map.drawmapscale(151.5, -34.73, 151, -34, 10, barstyle='simple', units='km',
                 fontsize=15, yoffset=None, labelstyle='simple', fontcolor='k',
                 fillcolor1='k', fillcolor2='k', ax=None, format='%d',
                 zorder=None)

plt.title("Sydney Metropolitan Region \n From MEGAN files", fontsize=20)
