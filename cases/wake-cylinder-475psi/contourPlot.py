from fluidfoam import readmesh 
from fluidfoam import readvector
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np

path = '/home/andy/OpenFOAM/andy-v1706/run/rarefiedPlume/cases/wake-cylinder'

# coordinates 
x, y =  readmesh(path, structured = False, boundary = 'inflow')

# Read velocity field
vel = readvector(path, '0.0495', 'U_Ar', False)

# Re-arange coordinates 
x = x.reshape(240, 920)
y = y.reshape(240, 920)
y = y[::-1, :]

# Magnitude of velocity field
magU = np.sqrt( vel[0, :]**2 + vel[1, :]**2)
magU.resize(240, 920)
magU = magU[::-1, :]

# Circle represents the turbine
circle = Circle((0,0), 5.0, color = 'black', fill = False)

# Plot the contours
fig = plt.figure(figsize = (6, 1.5))
plt.rc('font', size  = 18)
plt.rcParams["font.family"] = "monospace"
plt.title('Test')
plt.xlabel('x')
plt.ylabel('y')

ax = fig.add_subplot(1, 1, 1)

img = ax.imshow(
    magU,
    interpolation = 'bilinear',
    cmap = 'jet',
    extent = [x.min(), x.max(), y.min(), y.max()]
)

cbar = fig.colorbar(
    ax = ax,
    mappable = img,
    orientation = 'horizontal',
    pad = 0.3,
    label = 'm/s'
)

ax.add_artist(circle)
plt.savefig(path + '/contour.png')