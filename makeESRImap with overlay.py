# -*- coding: utf-8 -*-
"""
Create map with ESRI API Imagery 

@author: Jordan Capnerhurst 2016
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from mpl_toolkits.basemap import maskoceans

# Set colour map
cool = cm = plt.get_cmap('jet')

f = netCDF4.Dataset('O:/Honours_data/OEHYear/CTM_LAI_SoilType D3/ccam_3000m_lai_apr_RasterToN.nc', 'r')
f1 = netCDF4.Dataset('O:/Honours_data/OEHYear/CTM_LAI_SoilType D3/soiltype_3000m_RasterToNetCD.nc', 'r')
# [ time, source, lon, lat ]
# plot Daily average
vlai2 = f.variables['ccam_3000m_lai_apr.grd'][15:60, 0:42]
vsoil2 = f1.variables['soiltype_3000m.grd'][15:55, 0:42]


datalats = f.variables['y'][15:60] 
datalons = f.variables['x'][0:42]

mlons, mlats = np.meshgrid(datalons, datalats)

plt.figure(figsize=(15, 15))

map = Basemap(llcrnrlon=150.1,
              llcrnrlat=-34.7177,
              urcrnrlon=151.651,
              urcrnrlat=-33.5651, epsg=4269)
              
    #http://server.arcgisonline.com/arcgis/rest/services
    #EPSG Number of America is 4269

# Mask oceans and dams
mocedata = maskoceans(mlons, mlats, vlai2,  inlands=True, resolution='h',
                      grid=1.25)

# Colour mesh
map.pcolormesh(mlons, mlats, mocedata, latlon=True, zorder=1, vmax= 6,
               cmap=cool, alpha=0.1)
               
# Add colour bar
col = map.pcolormesh(mlons, mlats, mocedata,
                     latlon=True, zorder=0.45, vmax= 6, 
                     cmap=cool, alpha=0.1)  # set limits on colour scales

cb = map.colorbar(col, "right", size="5%", pad="2%")
cb.set_label('Average Total Biogenic Emissions (kg/$km^2$/hr)', fontsize=10)

plt.title("CTM Domain 3 3x3 Average Monthly Distribution of Biogenic Emissions During Daylight Hours\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)

map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)
plt.show()