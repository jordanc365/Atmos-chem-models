# -*- coding: utf-8 -*-
"""
Create total emissions map CTM 2013 

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

# Import data
m1 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JAN.nc', 'r')  # January
m2 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/FEB.nc', 'r')  # Febuary
m3 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/MAR.nc', 'r')  # March
m4 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/APR.nc', 'r')  # April
m5 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/MAY.nc', 'r')  # May
m6 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JUN.nc', 'r')  # June
m7 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JUL.nc', 'r')  # July
m8 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/AUG.nc', 'r')  # August
m9 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/SEP.nc', 'r')  # September
m10 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/OCT.nc', 'r')  # October
mnov = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/NOV.nc', 'r')  # November
mdec = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/DEC.nc', 'r')  # Decmber


# Import variables
# Store_bio

# [ time, source, lon, lat ] Jan
# plot Daily average
m11 = m1.variables['store_Bio'][:, 14:24, 0, :, :]
m12 = m1.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jan = np.concatenate((m11, m12), axis=1)/81

# average
janav = np.mean(jan, axis=(0, 1))
jantot = np.sum(jan, axis=(0, 1))
print jantot/1000
# Reshape for month
janres = jan.reshape(744, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] Feb
# plot Daily average
m21 = m2.variables['store_Bio'][:, 14:24, 0, :, :]
m22 = m2.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
feb = np.concatenate((m21, m22), axis=1)/81
febtot = np.sum(feb, axis=(0, 1,2,3))
print febtot/1000
# average
febav = np.mean(feb, axis=(0, 1))
febtot = np.sum(feb, axis=(0, 1))
# Reshape for month
febres = feb.reshape(672, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] Mar
# plot Daily average
m31 = m3.variables['store_Bio'][:, 14:24, 0, :, :]
m32 = m3.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
mar = np.concatenate((m31, m32), axis=1)/81

# average
marav = np.mean(mar, axis=(0, 1))
martot = np.sum(mar, axis=(0, 1))
print martot/1000

# Reshape for month
marres = mar.reshape(744, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] Apr
# plot Daily average
m41 = m4.variables['store_Bio'][:, 14:24, 0, :, :]
m42 = m4.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
apr = np.concatenate((m41, m42), axis=1)/81

# average
aprav = np.mean(apr, axis=(0, 1))
aprtot = np.sum(apr, axis=(0, 1))
print aprtot/1000

# Reshape for month
aprres = apr.reshape(720, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] May
# plot Daily average
m51 = m5.variables['store_Bio'][:, 14:24, 0, :, :]
m52 = m5.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
may = np.concatenate((m51, m52), axis=1)/81

# average
mayav = np.mean(may, axis=(0, 1))
maytot = np.sum(may, axis=(0, 1))
print maytot/1000

# Reshape for month
mayres = may.reshape(744, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] june
# plot Daily average
m61 = m6.variables['store_Bio'][:, 14:24, 0, :, :]
m62 = m6.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jun = np.concatenate((m61, m62), axis=1)/81

# average
junav = np.mean(jun, axis=(0, 1))
juntot = np.sum(jun, axis=(0, 1))
print juntot/1000

# Reshape for month
junres = jun.reshape(720, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] July
# plot Daily average
m71 = m7.variables['store_Bio'][:, 14:24, 0, :, :]
m72 = m7.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jul = np.concatenate((m71, m72), axis=1)/81

# average
julav = np.mean(jul, axis=(0, 1))
jultot = np.sum(jul, axis=(0, 1))
print jultot/1000

# Reshape for month
julres = jul.reshape(720, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] aug
# plot Daily average
m81 = m8.variables['store_Bio'][:, 14:24, 0, :, :]
m82 = m8.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
aug = np.concatenate((m81, m82), axis=1)/81

# average
augav = np.mean(aug, axis=(0, 1))
augtot = np.sum(aug, axis=(0, 1))
print augtot/1000

# Reshape for month
augres = aug.reshape(744, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] sep
# plot Daily average
m91 = m9.variables['store_Bio'][:, 14:24, 0, :, :]
m92 = m9.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
sep = np.concatenate((m91, m92), axis=1)/81

# average
sepav = np.mean(sep, axis=(0, 1))
septot = np.sum(sep, axis=(0, 1))
print septot/1000

# Reshape for month
sepres = sep.reshape(720, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] oct
# plot Daily average
m101 = m10.variables['store_Bio'][:, 14:24, 0, :, :]
m102 = m10.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
oct = np.concatenate((m101, m102), axis=1)/81

# average
octav = np.mean(oct, axis=(0, 1))
octtot = np.sum(oct, axis=(0, 1))
print octtot/1000

# Reshape for month
octres = oct.reshape(696, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] nov
# plot Daily average
m111 = mnov.variables['store_Bio'][:, 14:24, 0, :, :]
m112 = mnov.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
nov = np.concatenate((m111, m112), axis=1)/81

# average
novav = np.mean(nov, axis=(0, 1))
novtot = np.sum(nov, axis=(0, 1))
print novtot/1000

# Reshape for month
novres = nov.reshape(720, 60, 60) 

##############################################################################

# [ time, source, lon, lat ] dec
# plot Daily average
m121 = mdec.variables['store_Bio'][:, 14:24, 0, :, :]
m122 = mdec.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
dec = np.concatenate((m121, m122), axis=1)/81

# average
decav = np.mean(dec, axis=(0, 1))
dectot = np.sum(dec, axis=(0, 1))
print dectot/1000

# Reshape for month
decres = dec.reshape(720, 60, 60) 

##############################################################################
# concatanate all months 

year = np.concatenate((janres,febres,marres,aprres,mayres,junres,julres,augres,sepres,octres,novres,decres))

totdf = [sum(x) for x in zip(jantot, febtot, martot, aprtot, maytot, juntot, jultot, augtot, septot, octtot, novtot, dectot)]
totarrdf = np.array(totdf)/1000  # milimoles to kilograms

##############################################################################
# Define map parameters

max = 37.864
min = 0
trans = 0.2

##############################################################################

# Other variables
at = (m2.variables['lndtype'][0, :, :])
st = (m2.variables['soiltype'][0, :, :])
lai = (m1.variables['lai'][0, :, :])
t = (m2.variables['temp_a'][:, 0:24, 0, :, :])
g = (m2.variables['skin_temp'][:, 0:24, :, :])


'''
# Mean over array
v1 = m2con.mean(axis=(0, 1))
t1 = t.mean(axis=(0, 1))
g1 = g.mean(axis=(0, 1))
#cc1 = cc.mean(axis=(0, 1))
#cr1 = cr.mean(axis=(0, 1))
#ft1 = cr.mean(axis=(0, 1))
'''

datalats = m2.variables['lat'][:] 
datalons = m2.variables['lon'][:]

# Map
#  for 'low', not a numeral 1
map = Basemap(llcrnrlon=147.804, llcrnrlat=-36.7246,  # Lower left corner
              urcrnrlon=153.114, urcrnrlat=-31.4146, epsg=4269)  # Upper Right corner


plt.figure(figsize=(10, 10))

mlons, mlats = np.meshgrid(datalons, datalats)

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, totarrdf,  inlands=True, resolution='h',
                      grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, latlon=True, zorder=1, vmin=min, vmax=max,
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
'''
# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 1), labels=[0, 0, 0, 1], fontsize=10,
                  color='black', linewidth=2)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=10,
                  color='black', linewidth=2)


# Add colour bar
col = map.pcolormesh(mlons, mlats, mocedata,
                     latlon=True, zorder=1, vmin=min, vmax=max, alpha=trans,
                     cmap=cool)  

cb = map.colorbar(col, "bottom", size="5%", pad="2%")
cb.set_label('Totoal BVOC Emissions (Tonnes/year/$km^2$)', fontsize=15)

#plt.title("CTM Domain 3 3x3 Leaf Area Index\
#\n Sydney Metropolitan Region Feburary 2011", fontsize=10)


map.arcgisimage(service='World_Topo_Map', xpixels = 400, verbose= True)

plt.show()
# plt.savefig('./bmap_syd.png')
