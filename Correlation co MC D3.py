# -*- coding: utf-8 -*-
"""
Calculate correlation CSIRO CTM MEGAN

@author: Jordan Capnerhurst 2016
"""

# import things

import netCDF4
import numpy as np
from scipy.stats.stats import pearsonr

f0 = netCDF4.Dataset('O:/Honours_data/KMEGAN/KFebmegv2.nc',
                    'r')

f = netCDF4.Dataset('O:/Honours_data/KMEGAN/SPS1_3km_pft.nc','r') 

flai = netCDF4.Dataset('O:/Honours_data/KMEGAN/syd_3km_laifeb.nc','r')    

ftemp = netCDF4.Dataset('O:/Honours_data\KMEGAN\FebMEG.nc',
                    'r')             

flai2 = netCDF4.Dataset('O:/Honours_data/OEHYear/CTM_LAI_SoilType D3/ccam_3000m_lai_apr_RasterToN.nc', 'r')

# import variables

v = f.variables['PFT'][4, 3:44, 1:43]

v2 = f.variables['PFT'][6, 3:44, 1:43]

sv = (f0.variables['store_Megan'][1:27, :, 0, 3:44, 1:43])  # Isoprene
sv2 = (f0.variables['store_Megan'][1:27, :, 1, 3:44, 1:43])  # Monoterp

vlai = flai.variables['lai'][3:44, 1:43]

t = ftemp.variables['temp_a'][0:28, 0:24, 0, 3:44, 1:43]*1.025
g = ftemp.variables['skin_temp'][0:28, 0:24, 3:44, 1:43]*1.025

vlai2 = flai2.variables['ccam_3000m_lai_apr.grd'][ 0:41, 0:42]


# Mean over array
vm = sv.mean(axis=(0, 1))
ter1 = sv2.mean(axis=(0, 1))

t1 = t.mean(axis=(0, 1))
g1 = g.mean(axis=(0, 1))

# Add broad leaf PFTs
totdf = [sum(x) for x in zip(v, v2)]
totarrdf = np.array(totdf)  # milimoles to kilograms

# Add monoterpenes and isoprene to same array daily
totdf2 = [sum(x) for x in zip(vm, ter1)]
totarrdf2 = np.array(totdf2)/9  # milimoles to kilograms


rv1 = np.reshape(vlai2, 1722)

rv2 = np.reshape(totarrdf2, 1722)

corr = pearsonr(rv1, rv2)

print'correlation coefficient:', corr