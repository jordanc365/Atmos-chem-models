# -*- coding: utf-8 -*-
"""
Calculate spatial correlation MEGAN offline

@author: Jordan Capnerhurst 2016
"""

# import 

import netCDF4
import numpy as np
from scipy.stats.stats import pearsonr

f0 = netCDF4.Dataset('O:/Honours_data/JMEGAN/FebMegD3.nc',
                    'r')
f = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/btr200121_30sec.nc','r') #broad leaf tree
f1 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/hrb200121_30sec.nc', 'r') #Herbacioius
f2 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/ntr200121_30sec.nc', 'r') # needle leaf
f3 = netCDF4.Dataset('O:/Honours_data/JMEGAN/Plants/shr200121_30sec.nc', 'r') # Shrubs

frad = netCDF4.Dataset('O:/Honours_data/JMEGAN/rad1.nc','r')

flai = netCDF4.Dataset('O:/Honours_data/JMEGAN/LAI/laiv200302_30sec.nc',
                    'r')

# import variables

v0 = (f0.variables['ISOP'][:, :, 0,  :, 3:151])*68.12 # Isoprene
ter = (f0.variables['TERP'][:, :, 0, :, 3:151])*136.298 # terpenes

v = (f.variables['Broadleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'][1837:1965,6012:6160])
v1 = (f1.variables['Herbaceous_vegetation_cover_fraction_for_year_2001_(m2_per_m2)'][1837:1965,6012:6160])
v2 = (f2.variables['Needleleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'][1837:1965,6012:6160])
v3 = (f3.variables['Shrub_cover_fraction_for_year_2001_(m2_per_m2)'][1837:1965,6012:6160])

# [ time, source, lon, lat ]
rad = (frad.variables['GSW'][:, 0:23, 0, :, 3:151]) # Radiation at ground level
temp = (frad.variables['TEMP2'][:, 0:23, 0, :, 3:151]) # temp at 2m

LAI = (flai.variables['LAI_for_Feb_2003_(m2_per_m2)'][1837:1965,6012:6160])*0.001 # LAI

# Mean over array
vm = v0.mean(axis=(0, 1))
ter1 = ter.mean(axis=(0, 1))

radav = rad.mean(axis=(0, 1))
tempav = temp.mean(axis=(0, 1))

# Add monoterpenes and isoprene to same array daily
totdf = [sum(x) for x in zip(vm, ter1)]
totarrdf = np.array(totdf)*3.6/9  # milimoles to kilograms


rv1 = np.reshape(totarrdf, 18944)

rv2 = np.reshape(LAI, 18944)

corr = pearsonr(st1, totarr)

print'correlation coefficient:', corr