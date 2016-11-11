# -*- coding: utf-8 -*-
"""
Create PFT map MEGAN

@author: Jordan Capnerhurst 2016
"""

# import things
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4

# Set colour map
cool = cm = plt.get_cmap('jet')

f = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/btr200121_30sec.nc','r') #broad leaf tree
f1 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/hrb200121_30sec.nc', 'r') #Herbacioius
f2 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/ntr200121_30sec.nc', 'r') # needle leaf
f3 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/shr200121_30sec.nc', 'r') # Shrubs

# [ time, source, lon, lat ]
v = (f.variables['Broadleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'][1831:1972,6012:6160])
vma = np.ma.masked_less_equal(v, 10)

v1 = (f1.variables['Herbaceous_vegetation_cover_fraction_for_year_2001_(m2_per_m2)'][1831:1972,6012:6160])
vma1 = np.ma.masked_less_equal(v1, 10)

v2 = (f2.variables['Needleleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'][1831:1972,6012:6160])
vma2 = np.ma.masked_less_equal(v2, 10)

v3 = (f3.variables['Shrub_cover_fraction_for_year_2001_(m2_per_m2)'][1831:1972,6012:6160])
vma3 = np.ma.masked_less_equal(v3, 10)

datalons = f.variables['lon'][6012:6160]
datalats = f.variables['lat'][1831:1972]

map = Basemap(llcrnrlon=150.1,
              llcrnrlat=-34.7177,
              urcrnrlon=151.651,
              urcrnrlat=-33.5651, epsg=4269)
              
    #http://server.arcgisonline.com/arcgis/rest/services
    #EPSG Number of America is 4269

plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# set transparancy
trans = 1

# Colour mesh
map.pcolormesh(mlons, mlats, vma3, latlon=True, zorder=1, vmax=100, alpha=trans)

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
col = map.pcolormesh(mlons, mlats, vma3, latlon=True, zorder=1, cmap=cool, vmax=100, alpha=trans)
cb = map.colorbar(col, "bottom", size="5%", pad="2%")
cb.set_label('Percentage Cover', fontsize=15)

plt.title("MEGAN Domain 3 Shrub cover \n Sydney Metropolitan \
Region Feburary 2011", fontsize=10)

map.arcgisimage(service='World_Imagery', xpixels = 1000, verbose= True)
plt.show()
#plt.savefig('./bmap_syd.tiff')


