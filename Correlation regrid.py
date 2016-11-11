# -*- coding: utf-8 -*-
"""
Regridded spatial representation of temporal correlation CTM 2011

@author: Jordan Capnerhurst 2016
"""

# import things
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from mpl_toolkits.basemap import maskoceans
import scipy.interpolate

# Set colour map
cool = cm = plt.get_cmap('jet')

#f = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc', 'r')
# [ time, source, lon, lat ]
# plot Daily average


f = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/FEB.nc','r')

# new 10-km grid
glon3 = f.variables['lon'][:]
glat3 = f.variables['lat'][:]

# dlon, dlat = spacing (keep extra just in case)
dlon3 = glon3[-1] - glon3[-2]
dlat3 = glat3[-1] - glat3[-2]


f = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc','r')

# old 1-km grid
glon1 = f.variables['lon'][7:67]
glat1 = f.variables['lat'][7:67]

v = f.variables['store_Bio'][:, :, 0, 7:67, 7:67]
v2 = f.variables['store_Bio'][:, 0:14, 0, 6:46, 6:48]

# concatenate back to 24 hour period
#vcon = np.concatenate((v, v2), axis=1)

# Select desired period 
#vhalf= vcon[:, :, :, :]
#vhalf2= vcon[:, 20:23, :, : ]

#vcon2 = np.concatenate((vhalf, vhalf2), axis=1)



at = (f.variables['lndtype'][9, 6:46, 6:48])
st = (f.variables['soiltype'][0, :, :])
t = f.variables['temp_a'][0:28, 0:24, 0,  7:67, 7:67]
g = f.variables['skin_temp'][0:28, 0:24, :, :]
cc = f.variables['cloud_cover'][0:28, 0:24, 0, :, :]
cr = f.variables['cwater'][0:28, 0:24, 0, :, :]
z = f.variables['z'][0, :, :]
ft = f.variables['thstr'][0:28, 0:24, :, :]
lai = (f.variables['lai'][0, 6:46, 6:48])

#Add emissions to get total emissions
#total = np.sum(v, axis=(0, 1, 2, 3))*0.91

# Mean over array
#v1 = vhalf.mean(axis=(0, 1))/9
t1 = t.mean(axis=(0))
g1 = g.mean(axis=(0, 1))
cc1 = cc.mean(axis=(0, 1))
cr1 = cr.mean(axis=(0, 1))
ft1 = cr.mean(axis=(0, 1))

'''
# Correlation using for loops temp and bio
dat = np.zeros((v.shape[2], v.shape[3]))
for ii in range(v.shape[2]):
        for jj in range(v.shape[3]):
            vtmp = v[:, :, ii, jj].mean(axis=(0, 1))
            ttmp = t[:, :, ii, jj].mean(axis=(0, 1))
            dat[ii, jj] = np.correlate([vtmp], [ttmp])/850.543823242
'''         
# Correlation using for loops temp and LAI
dat2 = np.zeros((v.shape[2], v.shape[3]))
for ii in range(v.shape[2]):
        for jj in range(v.shape[3]):
            vtmp = v[:, :, ii, jj].mean(axis=(0, 1))
            ttmp = vflip[ii, jj]#.mean(axis=(0, 1))
            dat2[ii, jj] = np.correlate([vtmp], [ttmp])
            
            datr = dat2




LAI = datr


# Cut LAI & old lon/lat to match size of new grid (makes it a lot faster!!)
lonind = np.where((glon1 >= glon3[0]-dlon3) & (glon1 <= glon3[-1]+dlon3))
LAI = LAI[:,lonind[0][:]]
glon1 = glon1[lonind[0][:]]

latind = np.where((glat1 >= glat3[0]-dlat3) & (glat1 <= glat3[-1]+dlat3))
LAI = LAI[latind[0][:],:]
glat1 = glat1[latind[0][:]]

#create mesh
X, Y = np.meshgrid(glon1, glat1)
XI, YI = np.meshgrid(glon3, glat3)

#interpolate
LAInew=scipy.interpolate.griddata((X.flatten(),Y.flatten()),LAI.flatten() , (XI,YI),method='linear')



datalats = f.variables['lat'][7:67] 
datalons = f.variables['lon'][7:67]

# Map
#  for 'low', not a numeral 1
map = Basemap(llcrnrlon=150.1,
              llcrnrlat=-34.7177,
              urcrnrlon=151.651,
              urcrnrlat=-33.5651, epsg=4269)
    #http://server.arcgisonline.com/arcgis/rest/services
    #EPSG Number of America is 4269

plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# set transparancy / max and min
trans = 0.2
max = 1
min = 0

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, LAInew,  inlands=True, resolution='h',
                      grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, latlon=True, zorder=1,
               cmap=cool, alpha=trans, vmax = max, vmin = min)
'''
#map.drawstates(color='black', linewidth=3)
#map.drawcoastlines(color='black', linewidth=3)
# map.drawcountries()
# map.fillcontinents('white')
map.drawmapboundary()
# map.drawrivers(color='black', linewidth=2)
# map.shadedrelief()


# City markers and names
#lons = [151.2070, 150.8931, 151.7789, 149.1287]
#lats = [-33.8675, -34.4250, -32.9267, -35.2820]
#x, y = map(lons, lats)
#map.plot(x, y, 'ro', markersize=15)

labels = ['   Sydney', '   Wollongong', '   ', '']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt+5000, ypt-9000, label, color='black', fontsize=20)
'''
#Draw Meridians and labels
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
col = map.pcolormesh(mlons, mlats, mocedata,
                     latlon=True, zorder=0.45, 
                     cmap=cool, alpha=trans, vmax = max, vmin = min)  # set limits on colour scales

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Correlation Between LAI and Emissions', fontsize=15)


plt.title("CTM Domain 3 3x3 Average Monthly Distribution of Biogenic Emissions During Daylight Hours\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)

map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)

plt.show()
# plt.savefig('./bmap_syd.png')
