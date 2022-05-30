# Author: Andy Torres
# Purpose: Load inflow data from .dat file and plot relevant information using matplotlib. 


import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np


# Declare input data, angles from 0 to 90 degrees
theta = np.linspace(0, 1.57, 33)
print(theta)

radius = 1.05
x = np.round(radius * np.cos(theta), 2)
print(x)
y = np.round(radius * np.sin(theta), 2)
print(y)

velocity = 2
u = np.round(velocity* np.cos(theta), 2)
print(u)

v = np.round(velocity* np.sin(theta), 2)
print(v)

# Graph Velcoity data as a verctor function, save to file. 
plt.quiver(x, y, u, v)
plt.legend()
curr_dir_path = os.path.dirname(os.path.realpath(__file__))
plt.savefig(curr_dir_path + '/inflowCoordinates.png')
plt.close()

# Graph Nd data as a function of (x,y)

