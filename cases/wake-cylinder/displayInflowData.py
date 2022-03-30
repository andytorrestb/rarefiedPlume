# Author: Andy Torres
# Purpose: Load inflow data from .dat file and plot relevant information using matplotlib. 


import matplotlib.pyplot as plt
import pandas as pd
import os

# Load data, assign columns
df = pd.read_table("inflowData.dat", sep = "\s+")
df.columns = ["x coordinate", "x velocity", "y velocity", "temperature", "Nd"]

# Add y - coordinates to data set.

# Equation of a circle
# (x-h)^2 + (y-k)^2 = r^2

# Solving for y = f(x) yields,
# y = k + [r^2 -(x-h)^2] ^ 0.5

k = 0
h = -0.37465 # 11.75 + 3 inches converted to meters, in negative direction
r = 0.1524 # 6 inches converted to meters
df["y coordinate"] = k + pow(r**2 -pow(df["x coordinate"]-h, 2), 0.5)

# Re- order columns
df = df[["x coordinate", "y coordinate", "x velocity", "y velocity", "temperature", "Nd"]]
print(df)

# Graph Velcoity data as a verctor function, save to file. 
plt.quiver(df["x coordinate"] , df["y coordinate"], df["x velocity"], df["x velocity"])
plt.legend()
curr_dir_path = os.path.dirname(os.path.realpath(__file__))
plt.savefig(curr_dir_path + '/inflowCoordinates.png')
plt.close()

# Graph Nd data as a function of (x,y)

