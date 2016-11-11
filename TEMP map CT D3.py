# -*- coding: utf-8 -*-
"""
Create temperature map CTM 2011

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

f = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc', 'r')
# [ time, source, lon, lat ]
# plot Daily average
t = f.variables['temp_a'][0:28, 0:23, 0, :, :]
g = f.variables['skin_temp'][0:28, 0:23, 6:46, 6:48]
cc = f.variables['cloud_cover'][0:28, 0:23, 0, 6:46, 6:48]
v = f.variables['store_Bio'][0:28, 0:23, 0, 6:46, 6:48]
#Add emissions to get total emissions
#total = np.sum(v, axis=(0, 1, 2, 3))*0.91

# Mean over array
t1 = t.mean(axis=(0 ,1))
v1 = v.mean(axis=(0 ,1))
g1 = g.mean(axis=(0, 1))
cc1 = cc.mean(axis=(0, 1))

'''
# Correlation using for loops
dat = np.zeros((v.shape[2], v.shape[3]))
for ii in range(v.shape[2]):
    for jj in range(v.shape[3]):
        vtmp = v[:, :, ii, jj].mean(axis=(0, 1))
        ttmp = t[:, :, ii, jj].mean(axis=(0, 1))
        dat[ii, jj] = np.correlate([vtmp], [ttmp])
'''
    
datalats = f.variables['lat'][:] 
datalons = f.variables['lon'][:]

map = Basemap(llcrnrlon=150.1,
              llcrnrlat=-34.7177,
              urcrnrlon=151.651,
              urcrnrlat=-33.5651, epsg=4269)

plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

trans= 0.2

# Mask oceans and dams
#mocedata = maskoceans(mlons, mlats, t1,  inlands=True, resolution='h',
                      #grid=1.25)
# Colour mesh
map.pcolormesh(mlons, mlats, t1, latlon=True, zorder=1,vmax=29, vmin=19,
               cmap=cool, alpha=trans)

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
col = map.pcolormesh(mlons, mlats, t1, vmax=29, vmin=19,
                     latlon=True, zorder=1, 
                     cmap=cool, alpha=trans)  # set limits on colour scales

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Ambient Temperature ($^o$c)', fontsize=15)

plt.title("CTM Domain 3 3x3 Leaf Area Index\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)

map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)
plt.show()
# plt.savefig('./bmap_syd.png')
