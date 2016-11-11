# -*- coding: utf-8 -*-
"""
Create map with multiple transparent layers 

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
v = (f.variables['store_Bio'][:, :, 0, :, :])/9.0
at = (f.variables['lndtype'][0, :, :])
st = (f.variables['soiltype'][0, :, :])
t = f.variables['temp_a'][0:28, 0:24, 0, :, :]
g = f.variables['skin_temp'][0:28, 0:24, :, :]
cc = f.variables['cloud_cover'][0:28, 0:24, 0, :, :]
cr = f.variables['cwater'][0:28, 0:24, 0, :, :]
z = f.variables['z'][0, :, :]
ft = f.variables['thstr'][0:28, 0:24, :, :]
lai = (f.variables['lai'][0, :, :])

# Add emissions to get total emissions
#total = np.sum(v, axis=(0, 1))*9*0.91

# Mean over array
v1 = v.mean(axis=(0, 1))*0.91
t1 = t.mean(axis=(0))
g1 = g.mean(axis=(0, 1))
cc1 = cc.mean(axis=(0, 1))
cr1 = cr.mean(axis=(0, 1))
ft1 = cr.mean(axis=(0, 1))

# Correlation using for loops
dat = np.zeros((v.shape[2], v.shape[3]))
for ii in range(v.shape[2]):
    for jj in range(v.shape[3]):
        vtmp = v[:, :, ii, jj].mean(axis=(0, 1))
        ttmp = lai[ii, jj]#.mean(axis=(0, 1))
        dat[ii, jj] = np.correlate([vtmp], [ttmp])

datalons = f.variables['lon'][:]
datalats = f.variables['lat'][:]

# Map
#  for 'low', not a numeral 1
map = Basemap(projection='merc', lat_0=-33, lon_0=151,
              resolution='h', area_thresh=0.1,
              llcrnrlon=150.1, llcrnrlat=-34.7177,  # Lower left corner
              urcrnrlon=151.651, urcrnrlat=-33.5651)  # Upper Right corner


plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, lai,  inlands=True, resolution='h',
                      grid=1.25)
                      
mocedata2 = maskoceans(mlons, mlats, v1,  inlands=True, resolution='h',
                      grid=1.25)


# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, vmin=0, vmax=5.4, latlon=True,
               cmap=cool)
               
map.pcolormesh(mlons, mlats, mocedata2, vmin=0, vmax=5.4, latlon=True, zorder=10,
               cmap=cool, alpha=0.2)

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
col = map.pcolormesh(mlons, mlats, mocedata, vmin=0, vmax=5.4,
                     latlon=True,
                     cmap=cool)  # set limits on colour scales
                     
col2 = map.pcolormesh(mlons, mlats, mocedata2, vmin=0, vmax=5.4,
                     latlon=True, zorder=10,
                     cmap=cool, alpha=0.2)  # set limits on colour scales
        

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('LAI ($m^2$/$m^2$)', fontsize=15)

cb = map.colorbar(col2, "bottom", size="5%", pad="2%")
cb.set_label('Biogenic emissions (kg/$km^2$/hour)', fontsize=15)

plt.title("CTM Domain 3 3x3  LAI and Total Biogenic Emissions\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)


plt.show()
# plt.savefig('./bmap_syd.png')
