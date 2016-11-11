# -*- coding: utf-8 -*-
"""
Basic correlation coeeficient script

@author: Jordan Capnerhurst 2016
"""

# Import things
from scipy.stats.stats import pearsonr
import netCDF4
import numpy as np
f = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc',
                    'r')
# [ time, source, lon, lat ]
v = (f.variables['store_Bio'][0, 0:23, 0, 30, 30])
t = f.variables['temp_a'][0, 0:23, 0, 30, 30]
g = f.variables['skin_temp'][0, 0:23, 30, 30]


p1 = pearsonr(v, t)
p2 = pearsonr(v, g)
p3 = pearsonr(t, g)
p4 = np.mean(v/89, axis=0)
p5 = np.mean(t, axis=0)
p6 = np.mean(g, axis=0)

print'Daily data 1/2/2013 at pt. 30,30 at 00 to 24 UTC'
print'correlation coefficient bio and ambient temp:', p1
print'correlation coefficient bio and skin temp:', p2
print'correlation coefficient ambient temp and skin temp:', p2
print'mean bio:', p4
print'mean ambient temp:', p5
print'mean skin temp:', p6

# data for monthly values
f2 = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc', 'r')

v2 = f2.variables['store_Bio'][0:27, 0:24, 0, :, :]
t2 = f2.variables['temp_a'][0:27, 0:24, :, :, :]
g2 = f2.variables['skin_temp'][0:27, 0:24, :, :]

# Average over lat and lon
v21 = v2.mean(axis=(2, 3))
t21 = t2.mean(axis=(2, 3, 4))
g21 = g2.mean(axis=(2, 3))


# Reshape array
r = np.reshape(v21, 648)
tr = np.reshape(t21, 648)
gr = np.reshape(g21, 648)

p11 = pearsonr(r, tr)
p21 = pearsonr(r, gr)
p31 = pearsonr(tr, gr)
p41 = np.mean(r/89, axis=0)
p51 = np.mean(tr, axis=0)
p61 = np.mean(gr, axis=0)

print'Monthly mean data'
print'correlation coefficient bio and ambient temp:', p11
print'correlation coefficient bio and skin temp:', p21
print'correlation coefficient ambient temp and skin temp:', p31
print'mean bio:', p41
print'mean ambient temp:', p51
print'mean skin temp:', p61

# Data for daily average values
f3 = netCDF4.Dataset('O:/Honours_data/KOEH/FebKOEH.nc', 'r')

# plot Daily average
v3 = f3.variables['store_Bio'][0:28, 0:24, 0, :, :]
t3 = f3.variables['temp_a'][0:28, 0:24, :, :]
g3 = f3.variables['skin_temp'][0:28, 0:24, :, :]

v31 = v3.mean(axis=(0, 2, 3))
t31 = t3.mean(axis=(0, 2, 3, 4))
g31 = g3.mean(axis=(0, 2, 3))


# Reshape array
r2 = v31
tr2 = t31
gr2 = g31

p13 = pearsonr(r2, tr2)
p23 = pearsonr(r2, gr2)
p33 = pearsonr(tr2, gr2)
p43 = np.mean(r2/89, axis=0)
p53 = np.mean(tr2, axis=0)
p63 = np.mean(gr2, axis=0)

print'Daily mean data'
print'correlation coefficient bio and ambient temp:', p13
print'correlation coefficient bio and skin temp:', p23
print'correlation coefficient ambient temp and skin temp:', p3
print'mean bio:', p43
print'mean ambient temp:', p53
print'mean skin temp:', p63
