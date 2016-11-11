# -*- coding: utf-8 -*-
"""
Create pcolourmesh CTM 2013 LAI and Soiltype

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

f = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_jan_RasterToN.nc', 'r')
f1 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_feb_RasterToN.nc', 'r')
f2 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_mar_RasterToN.nc', 'r')
f3 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_apr_RasterToN.nc', 'r')
f4 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_may_RasterToN.nc', 'r')
f5 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_jun_RasterToN.nc', 'r')
f6 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_jul_RasterToN.nc', 'r')
f7 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_aug_RasterToN.nc', 'r')
f8 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_sep_RasterToN.nc', 'r')
f9 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_oct_RasterToN.nc', 'r')
f10 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_nov_RasterToN.nc', 'r')
f11 = netCDF4.Dataset('O:/Honours_data/OEHYear/LAI D2 nc data/ccam_9000m_lai_dec_RasterToN.nc', 'r')



#f1 = netCDF4.Dataset('O:/Honours_data/OEHYear/CTM_LAI_SoilType D2/soiltype_3000m_RasterToNetCD.nc', 'r')
# [ time, source, lon, lat ]
# plot Daily average
v = f.variables['ccam_9000m_lai_jan.grd'][:, :]
vav = np.mean(v)
v1 = f1.variables['ccam_9000m_lai_feb.grd'][:, :]
vav1 = np.mean(v1)
v2 = f2.variables['ccam_9000m_lai_mar.grd'][:, :]
vav2 = np.mean(v2)
v3 = f3.variables['ccam_9000m_lai_apr.grd'][:, :]
vav3 = np.mean(v3)
v4 = f4.variables['ccam_9000m_lai_may.grd'][:, :]
vav4 = np.mean(v4)
v5 = f5.variables['ccam_9000m_lai_jun.grd'][:, :]
vav5 = np.mean(v5)
v6 = f6.variables['ccam_9000m_lai_jul.grd'][:, :]
vav6 = np.mean(v6)
v7 = f7.variables['ccam_9000m_lai_aug.grd'][:, :]
vav7 = np.mean(v7)
v8 = f8.variables['ccam_9000m_lai_sep.grd'][:, :]
vav8 = np.mean(v8)
v9 = f9.variables['ccam_9000m_lai_oct.grd'][:, :]
vav9 = np.mean(v9)
v10 = f10.variables['ccam_9000m_lai_nov.grd'][:, :]
vav10 = np.mean(v10)
v11 = f11.variables['ccam_9000m_lai_dec.grd'][:, :]
vav11 = np.mean(v11)

############################################################################
# plot average LAI over time
summervavyr =np.array([vav, vav1, vav11])
autumnvavyr =np.array([vav2,vav3, vav4])
wintervavyr =np.array([vav5, vav6, vav7])
springvavyr =np.array([vav8, vav9, vav10])
#months = np.arange(0, 12, 1)

#plt.plot(months, vavyr, linestyle='-', linewidth=2, c='r',
        # label='Average Ambient Temperature ($^o$c)')
         
plt.ylim(1, 3)
plt.xlim(0, 11)

# tick
plt.xticks(range(0, 12, 1 ), [str(i) for i in range(0, 12, 1)])
plt.yticks(np.arange(1, 3.1, 0.25))  # use floats for y value


############################################################################
#v1 = f1.variables['soiltype_3000m.grd'][:, :]

vflip3= np.flipud(v2)
vres =np.reshape(vflip, 3600)



datalats = f.variables['y'][:] 
datalons = f.variables['x'][:]

# Map
#  for 'low', not a numeral 1
map = Basemap(projection='merc', lat_0=-33, lon_0=151,
              resolution='h', area_thresh=0.1,
              llcrnrlon=147.804, llcrnrlat=-36.7246,  # Lower left corner
              urcrnrlon=153.114, urcrnrlat=-31.4146)  # Upper Right corner



plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, v4,  inlands=True, resolution='h',
                      grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, latlon=True, zorder=1,
               cmap=cool, vmax=8)
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
map.plot(x, y, 'ro', markersize=15)

labels = ['      Sydney', 'Wollongong', '  Newcastle', '   Canberra']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-90000, ypt+10000, label, color='white')
fontsize = 14

# North Arrow
lons = [152.7]
lats = [-36.5]
x, y = map(lons, lats)
map.plot(x, y, 'k^', markersize=30)

# Add scale bar
map.drawmapscale(151.5, -36.4, 151, -34, 50, barstyle='simple', units='km',
                 fontsize=15, yoffset=None, labelstyle='simple', fontcolor='k',
                 fillcolor1='k', fillcolor2='k', ax=None, format='%d',
                 zorder=None)

labels = ['North', '', '  ', '   ']
for label, xpt, ypt in zip(labels, x, y):
    plt.text(xpt-20000, ypt+30000, label)
fontsize = 40

# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 1), labels=[0, 0, 0, 1], fontsize=10,
                  color='black', linewidth=2)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=10,
                  color='black', linewidth=2)
# Add colour bar
col = map.pcolormesh(mlons, mlats, mocedata,
                     latlon=True, zorder=0.45, vmax=8, 
                     cmap=cool)  # set limits on colour scales

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Average Total Biogenic Emissions (kg/$km^2$/hr)', fontsize=10)

plt.title("CTM Domain 3 3x3 Average Monthly Distribution of Biogenic Emissions During Daylight Hours\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)


plt.show()
# plt.savefig('./bmap_syd.png')
