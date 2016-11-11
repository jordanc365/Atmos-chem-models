# -*- coding: utf-8 -*-
"""

Create plot from imported NETCDF data 

@author: Jordan Capnerhurst 2016
"""

import matplotlib.pyplot as plt

plt.plot()

# plot variables
f = "dataset"
plt.plot(f.variables['store_Bio'][1:28, 1, 0, 33, 33])

plt.yticks(np.arange(min(bios1), max(bios1)+1, 0.45))  # use floats for y value