# -*- coding: utf-8 -*-
"""
Create PFT map MEGAN

@author: Jordan Capnerhurst 2016
"""

# import things
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pyplot
import numpy as np
import netCDF4
import scipy.interpolate

# Set colour map
cool = cm = pyplot.get_cmap('rainbow')

f = netCDF4.Dataset('O:\Honours_data\KMEGAN\FebMEG.nc','r')

pf1 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/btr200121_30sec.nc','r') #broad leaf tree
pf2 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/hrb200121_30sec.nc', 'r') #Herbacioius
pf3 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/ntr200121_30sec.nc', 'r') # needle leaf
pf4 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/shr200121_30sec.nc', 'r') # Shrubs

# [ time, source, lon, lat ]
v = (pf1.variables['Broadleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'])
#vma = np.ma.masked_less_equal(v, 10)

v1 = (pf2.variables['Herbaceous_vegetation_cover_fraction_for_year_2001_(m2_per_m2)'])
#vma1 = np.ma.masked_less_equal(v1, 10)

v2 = (pf3.variables['Needleleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'])
#vma2 = np.ma.masked_less_equal(v2, 10)

v3 = (pf4.variables['Shrub_cover_fraction_for_year_2001_(m2_per_m2)'])
#vma3 = np.ma.masked_less_equal(v3, 10)

#------------------------------------------------------------------------------

# new 3-km grid
glon3 = f.variables['lon'][:]
glat3 = f.variables['lat'][:]

# dlon, dlat = spacing (keep extra just in case)
dlon3 = glon3[-1] - glon3[-2]
dlat3 = glat3[-1] - glat3[-2]


# old 1-km grid
glon1 = pf1.variables['lon'][:]
glat1 = pf1.variables['lat'][:]


# Cut LAI & old lon/lat to match size of new grid (makes it a lot faster!!)
lonind = np.where((glon1 >= glon3[0]-dlon3) & (glon1 <= glon3[-1]+dlon3))
v = v[:,lonind[0][:]]
glon1 = glon1[lonind[0][:]]

latind = np.where((glat1 >= glat3[0]-dlat3) & (glat1 <= glat3[-1]+dlat3))
v = v[latind[0][:],:]
glat1 = glat1[latind[0][:]]

#create mesh
X, Y = np.meshgrid(glon1, glat1)
XI, YI = np.meshgrid(glon3, glat3)

#interpolate
PFTnew=scipy.interpolate.griddata((X.flatten(),Y.flatten()),v.flatten() , (XI,YI), method='linear')

PFTmask = np.ma.masked_less_equal(PFTnew, 10)

#------------------------------------------------------------------------------
trans = 0.6

map = Basemap(llcrnrlon=150.1,
              llcrnrlat=-34.7177,
              urcrnrlon=151.651,
              urcrnrlat=-33.5651, epsg=4269)
              
    #http://server.arcgisonline.com/arcgis/rest/services
    #EPSG Number of America is 4269



f=pyplot.figure()

pyplot.figure(figsize=(10, 10))


    
# plot the data
col=map.pcolormesh(glon3,glat3,PFTmask,latlon=True,vmin=0, vmax = 100, alpha=trans,cmap=cool)


# Draw Meridians and labels
map.drawmeridians(np.arange(0, 360, 1), labels=[0, 0, 0, 1], fontsize=10,
                  color='black', linewidth=2)
map.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=10,
                  color='black', linewidth=2)

map.arcgisimage(service='World_Imagery', xpixels = 1000, verbose= True)

# add a color bar
cb = map.colorbar(col, "bottom", alpha=trans)

pyplot.title("MEGAN Domain 3 Shrub cover \n Sydney Metropolitan \
Region Feburary 2011", fontsize=10)



#plt.savefig('./bmap_syd.tiff')


