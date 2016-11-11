# -*- coding: utf-8 -*-
"""
Create pcolourmesh CTM LAI

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

f = netCDF4.Dataset('D:/2013 CCAM-CTM output/jan13/gmr_20130101_cb05_aer2_0.090dg_jan13.nc', 'r')

lai = (f.variables['lai'][:, :])

#Add emissions to get total emissions
#total = np.sum(v, axis=(0, 1, 2, 3))*0.91

# Mean over array
v1 = vhalf.mean(axis=(0, 1))*0.91/9
t1 = t.mean(axis=(0))
g1 = g.mean(axis=(0, 1))
cc1 = cc.mean(axis=(0, 1))
cr1 = cr.mean(axis=(0, 1))
ft1 = cr.mean(axis=(0, 1))

datalats = f.variables['lat'][:] 
datalons = f.variables['lon'][:]

# Map
#  for 'low', not a numeral 1
map = Basemap(projection='merc', lat_0=-33, lon_0=151,
              resolution='h', area_thresh=0.1,
              llcrnrlon=147.804, llcrnrlat=-36.7246,  # Lower left corner
              urcrnrlon=153.114, urcrnrlat=-31.4146)  # Upper Right corner


plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, lai,  inlands=True, resolution='h',
                      grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, latlon=True, zorder=1, 
               cmap=cool)
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
map.plot(x, y, 'ko', markersize=15)

labels = ['   Sydney', '   Wollongong', ' Newcastle', '    Canberra']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt+5000, ypt-9000, label, color='black', fontsize=20)

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 3), labels=[0, 0, 0, 1], fontsize=10)
map.drawparallels(np.arange(-90, 90, 3), labels=[1, 1, 1, 1], fontsize=10)

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

# Add colour bar
col = map.pcolormesh(mlons, mlats, mocedata,
                     latlon=True, zorder=1, 
                     cmap=cool)  # set limits on colour scales

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('LAI ($m^2$/$m^2$)', fontsize=15)

plt.title("CTM Domain 3 3x3 Leaf Area Index\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)


plt.show()
# plt.savefig('./bmap_syd.png')
