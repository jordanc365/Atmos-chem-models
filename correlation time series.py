# -*- coding: utf-8 -*-
"""
plot correlation over time 2013

@author: Jordan Capnerhurst
"""

import matplotlib.pyplot as plt
import numpy as np

temp = np.array([0.51, 0.82, 0.75, 0.75, 0.61, 0.44, 0.58, 0.49, 0.56, 0.67, 0.69, 0.73])

lai = np.array([0.98, 0.98, 0.98, 0.96, 0.95, 0.92, 0.91, 0.90, 0.90, 0.93, 0.95, 0.98])

month = np.arange(1, 13, 1)

plt.figure(figsize=(15, 5))

plt.plot(month, temp, linestyle='-', linewidth=2.0, c='r',
         label='Temperature')
         
plt.plot(month, lai, linestyle='-', linewidth=2.0, c='b',
         label='LAI')
         
plt.ylim(0, 1)
plt.xlim(1, 12)

plt.xticks(range(1, 13, 1 ), [str(i) for i in range(1, 13, 1)])
plt.yticks(np.arange(min(lai)-0.9, max(lai)+0.1, 0.1))  # use floats for y value

# Legend
plt.legend(loc='lower right', fontsize=10)