# -*- coding: utf-8 -*-
"""
Create map using ESRI API data 

@author: Jordan Capnerhurst 2016
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(10, 10))

map = Basemap(llcrnrlon=147.804, llcrnrlat=-36.7246,  # Lower left corner
              urcrnrlon=153.114, urcrnrlat=-31.4146, epsg=4269)
              
    # all map types can be found at http://server.arcgisonline.com/arcgis/rest/services
    #EPSG Number of America is 4269
'''
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
'''


# Add scale bar
#map.drawmapscale(151.5, -34.65, 151, -34, 5, barstyle='simple', units='km',
                 #fontsize=15, yoffset=None, labelstyle='simple', fontcolor='k',
                 #fillcolor1='k', fillcolor2='k', ax=None, format='%d', zorder=None)

map.arcgisimage(service='World_Topo_Map', xpixels = 706, verbose= True)

#labels = ['North', '', '  ', '   ']
#for label, xpt, ypt in zip(labels, x, y):
    #plt.text(xpt-10000, ypt+8000, label, color='black', fontsize=20)
    
plt.show()