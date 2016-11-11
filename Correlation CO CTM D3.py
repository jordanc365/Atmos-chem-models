# -*- coding: utf-8 -*-
"""
Calculate correlation Co. for any dataset 

@author: Jordan Capnerhurst 2016
"""

# import things

import netCDF4
import numpy as np
from scipy.stats.mstats import pearsonr

f = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc', 'r')

flai2 = netCDF4.Dataset('O:/Honours_data/OEHYear/CTM_LAI_SoilType D3/ccam_3000m_lai_apr_RasterToN.nc', 'r')


# import variables

v = f.variables['store_Bio'][:, :, 0, 6:46, 6:48]
at = (f.variables['lndtype'][5, 6:46, 6:48])
lai = (f.variables['lai'][0, 6:46, 6:48])
st = (f.variables['soiltype'][0, 6:46, 6:48])
t = f.variables['temp_a'][0:28, 0:24, :, 6:46, 6:48]
g = f.variables['skin_temp'][0:28, 0:24, 6:46, 6:48]
cc = f.variables['cloud_cover'][0:28, 0:24, 0, 6:46, 6:48]
z = f.variables['z'][0, 6:46, 6:48]

vlai2 = flai2.variables['ccam_3000m_lai_apr.grd'][ 15:55, 0:42]


# Mean over array
vm = v.mean(axis=(0, 1))


t1 = t.mean(axis=(0, 1))
g1 = g.mean(axis=(0, 1))
cc1 = cc.mean(axis=(0, 1))




rv1 = np.reshape(vm, 1680)

rv2 = np.reshape(vlai2, 1680)

corr = pearsonr(febtavres, febavres)

print'correlation coefficient:', corr