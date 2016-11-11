# -*- coding: utf-8 -*-
"""
Basic Numpy functions

@author: Jordan Capnerhurst 2016
"""

# import netCDF


import numpy as np


# Useful functions
r = 'imported_array'
r.size   # Size of array r
r.shape  # size of array r
ar = np.array(r)  # make data in r an array

np.amax(r)   # max of data r
np.amin(r)   # min of data r
