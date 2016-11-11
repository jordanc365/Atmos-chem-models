# -*- coding: utf-8 -*-
"""
Create 3x3km temperature map

@author: Jordan Capnerhurst 2016
"""

import numpy as np
import netCDF4
import scipy.interpolate
import matplotlib.pyplot as pyplot
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans

f = netCDF4.Dataset('O:\Honours_data\KMEGAN\FebMEG.nc','r')

# new 3-km grid
glon3 = f.variables['lon'][:]
glat3 = f.variables['lat'][:]

# dlon, dlat = spacing (keep extra just in case)
dlon3 = glon3[-1] - glon3[-2]
dlat3 = glat3[-1] - glat3[-2]



f = netCDF4.Dataset('O:/Honours_data/JMEGAN\Grid\GRIDCRO2D_d03.nc','r')
                      


# old 1-km grid
glon1 = f.variables['LON'][0, 0, 0, 10:155]
glat1 = f.variables['LAT'][0, 0, :, 0]

#------------------------------------------------------------------------------
f = netCDF4.Dataset('O:/Honours_data/JMEGAN/rad1.nc','r')
temp = (f.variables['TEMP2'][:, 14:24, 0, :, :]) # temp at 2m
temp2 = (f.variables['TEMP2'][:, 0:14, 0, :, :])

# concatenate back to 24 hour period
tercon = np.concatenate((temp, temp2), axis=1)


# Averages
TEMP = tercon.mean(axis=(0 ,1))-273.15

#------------------------------------------------------------------------------

# Cut LAI & old lon/lat to match size of new grid (makes it a lot faster!!)
lonind = np.where((glon1 >= glon3[0]-dlon3) & (glon1 <= glon3[-1]+dlon3))
TEMP = TEMP[:,lonind[0][:]]
glon1 = glon1[lonind[0][:]]

latind = np.where((glat1 >= glat3[0]-dlat3) & (glat1 <= glat3[-1]+dlat3))
TEMP = TEMP[latind[0][:],:]
glat1 = glat1[latind[0][:]]

#create mesh
X, Y = np.meshgrid(glon1, glat1)
XI, YI = np.meshgrid(glon3, glat3)

#interpolate
TEMPnew=scipy.interpolate.griddata((X.flatten(),Y.flatten()),TEMP.flatten() , (XI,YI),method='linear')

#mocedata = maskoceans(XI, YI, TEMPnew,  inlands=True, resolution='h',
                      #grid=1.25)

trans = 0.3




# Plot new
#----------

# Set up the map
# Set up the map
map = Basemap(llcrnrlon=150.1,
              llcrnrlat=-34.7177,
              urcrnrlon=151.651,
              urcrnrlat=-33.5651, epsg=4269 )

f=pyplot.figure()


pyplot.figure(figsize=(10, 10))

    
# plot the data
col=map.pcolormesh(glon3,glat3,TEMPnew,latlon=True, vmin=17.6, vmax = 24.8, alpha=trans)


# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 1), labels=[0, 0, 0, 1], fontsize=10,
                  color='black', linewidth=2)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=10,
                  color='black', linewidth=2)

map.arcgisimage(service='World_Topo_Map', xpixels = 820, verbose= True)


# add a color bar
cb = map.colorbar(col, "right", alpha=trans)
