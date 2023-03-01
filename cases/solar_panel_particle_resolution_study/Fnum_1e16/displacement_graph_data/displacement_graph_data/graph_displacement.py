import pandas as pd
import matplotlib.pyplot as plt

h20 = pd.read_csv('plume_pressure_distribution.csv')
h10 = pd.read_csv('plume_pressure_distribution_h10m.csv')
h05 = pd.read_csv('plume_pressure_distribution_h05m.csv')


# Label graph with bold characters
font_axis_publish = {
      'color':  'black',
      'weight': 'bold',
      'size': 22,
      }

# Graph Results
plt.title('Displacement Along the the Solar Panel (mm)', fontdict = font_axis_publish)
plt.ylabel('Z Displacement (mm)', fontdict = font_axis_publish)
plt.xlabel('X Displacement (mm)', fontdict = font_axis_publish)

plt.plot(h20['#Displacement(mm)'], h20['Displacement(mm)_Real'], label = 'Plune Impingement (h = 20 m)')
plt.plot(h10['#Displacement(mm)'], h10['Displacement(mm)_Real'], label = 'Plune Impingement (h = 10 m)')
plt.plot(h05['#Displacement(mm)'], h05['Displacement(mm)_Real'], label = 'Plune Impingement (h = 5 m)')
plt.legend()
plt.show()