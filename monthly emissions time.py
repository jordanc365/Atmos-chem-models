# -*- coding: utf-8 -*-
"""
Create monthly emissions timeseries 2013

@author: Jordan Capnerhurst 2016
"""

import numpy as np
from pylab import rcParams

year = np.array([5337.469, 4187.121, 3893.942, 1973.351, 856.073, 362.464, 
                 366.910, 626.229, 2132.338, 3249.037, 3942.992, 4355.3465])*9*0.000001

month = np.arange(0, 12, 1)

tot = np.sum(year)


rcParams['figure.figsize'] = 15, 5

plt.plot(month, year, linestyle='-', linewidth=2.0, c='b',
         label='Biogenic Emissions 2011')
         
# asthetics
plt.xlabel('Month')
plt.ylabel('Total Biogenic Emissions (Tg/ month)')
#plt.title('Monthly Ambiant temperature Feburary 2011 from CTM')
#plt.ylim(14, 26)
plt.xlim(0, 11)


# tick
plt.xticks(range(0, 12, 1), [str(i) for i in range(0, 12, 1)])
plt.yticks(np.arange(0, 0.051, 0.005))  # use floats for y value