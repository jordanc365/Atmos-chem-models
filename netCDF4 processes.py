# -*- coding: utf-8 -*-
"""
Usefull netCDF4 functions

@author: Jordan Capnerhurst 2016
"""

# import netCDF

import netCDF4
f = netCDF4.Dataset('O:/Honours_data/FebAVSA.nc', 'r')

# get list of variables
for v in f.variables:
    print(v)

# Acess a variable
f.variables['store_Bio'][:]

# Find dimensions of variable
f.variables['store_Bio'].dimensions

# Access variable in dimensions
f.variables['store_Bio'][1:24, 0, 0]

# Access variable in dimensions with date
f.variables['store_Bio'][1, 24, 0, 33, 59]
