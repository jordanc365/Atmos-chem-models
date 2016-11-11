# -*- coding: utf-8 -*-
"""
Create Emissions map MEGAN coupled to CTM

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

# Leaf area index used in model
flai = netCDF4.Dataset('O:/Honours_data/KMEGAN/syd_3km_laifeb.nc', 'r')

f = netCDF4.Dataset('O:/Honours_data/KMEGAN/KFebmegv2.nc',
                    'r')

at = (f.variables['lndtype'][0, :, :])
st = (f.variables['soiltype'][0, :, :])
t = f.variables['temp_a'][0:28, 0:24, 0, :, :]
g = f.variables['skin_temp'][0:28, 0:24, :, :]
cc = f.variables['cloud_cover'][0:28, 0:24, 0, :, :]
cr = f.variables['cwater'][0:28, 0:24, 0, :, :]
z = f.variables['z'][0, :, :]
ft = f.variables['thstr'][0:28, 0:24, :, :]
# lai = (f.variables['lai'][0, :, :])

# plot Daily average
sv = (f.variables['store_Megan'][1:27, 14:24, 0, :, :])  # Isoprene
sv2 = (f.variables['store_Megan'][1:27, 0:14, 0, :, :])  # Isoprene
sm = (f.variables['store_Megan'][1:27, 14:24, 1, :, :])  # Monoterp
sm2 = (f.variables['store_Megan'][1:27, 0:14, 1, :, :])  # Monoterp

# concatenate back to 24 hour period
vcon = np.concatenate((sv, sv2), axis=1)
tercon = np.concatenate((sm, sm2), axis=1)

vhalf= vcon[:, :, :, :]
#vhalf1= vcon[:, 20:23, :, :]

vhalf2= tercon[:, :, :, :]
#vhalf21= tercon[:, 20:23, :, :]

#vcon2 = np.concatenate((vhalf, vhalf1), axis=1)
#tercon2 = np.concatenate((vhalf2, vhalf21), axis=1)


lai = flai.variables['lai']  # leaf area index used model
lait = np.array(lai)

# Add emissions to get total emissions
#total = np.sum(v, axis=(1, 2))
#total1 = np.sum(ter, axis=(1, 2))

# Mean over array
v1 = vhalf.mean(axis=(0, 1))/9
ter1 = vhalf2.mean(axis=(0, 1))/9
t1 = t.mean(axis=(0, 1))
g1 = g.mean(axis=(0, 1))
cc1 = cc.mean(axis=(0, 1))*100
cr1 = cr.mean(axis=(0, 1))
ft1 = cr.mean(axis=(0, 1))

# Add monoterpenes and isoprene to same array 
totd = [sum(x) for x in zip(v1, ter1)]
totarrd = np.array(totd)/9 # milimoles to kilograms

megres =np.reshape(totarrd, 3600)

'''
# Correlation using for loops
dat = np.zeros((v.shape[2], v.shape[3]))
for ii in range(v.shape[2]):
    for jj in range(v.shape[3]):
        vtmp = totarrd[:, :, ii, jj].mean(axis=(0, 1))
        ttmp = lait[ii, jj]#.mean(axis=(0, 1))
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
max = 100
min = 0

# Mask oceans and dams
#mocedata = maskoceans(mlons, mlats, cc1, inlands=True, resolution='h',
                      #grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, cc1,  latlon=True, zorder=1,
               cmap=cool, alpha=trans, vmax=max, vmin=min)
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
col = map.pcolormesh(mlons, mlats, cc1, latlon=True, alpha=trans, vmax=max, vmin=min, cmap=cool)
'''
# Add scale bar
map.drawmapscale(151.5, -34.65, 151, -34, 5, barstyle='simple', units='km',
                 fontsize=15, yoffset=None, labelstyle='simple', fontcolor='k',
                 fillcolor1='k', fillcolor2='k', ax=None, format='%d',
                 zorder=None)
'''
cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Average Monoterpene Emissions (kg/$km^2$/hr)', fontsize=15)

plt.title("CTM Inline MEGAN Domain 3 3x3 Monthly Distribution of \
Biogenic Emissions \n Sydney Metropolitan Region  February 2011", fontsize=10)

map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)

plt.show()
# plt.savefig('./bmap_syd.png')
