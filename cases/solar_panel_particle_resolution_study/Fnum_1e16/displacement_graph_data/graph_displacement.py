import pandas as pd
import matplotlib.pyplot as plt

plume = pd.read_csv('plume_pressure_distribution.csv')
max_ = pd.read_csv('uniform_pressure_distribution_max.csv')
mean_ = pd.read_csv('uniform_pressure_distribution_avg.csv')
min_ = pd.read_csv('uniform_pressure_distribution_min.csv')

print(plume.dtypes)
print(plume.columns)
print(plume['#X'])
# plume.plot(kind = 'line', x = '#X', y='Displacement(mm)_Real')
# max_.plot(kind = 'scatter', x = '#X', y='Displacement(mm)_Real')
# mean_.plot(kind = 'scatter', x = '#X', y='Displacement(mm)_Real')
# min_.plot(kind = 'scatter', x = '#X', y='Displacement(mm)_Real')
# plt.show()
print(max_['#X'])
print(mean_['#X'])
print(min_['#X'])
print()
print(type(plume))

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

plt.plot(plume['#X'], plume['Displacement(mm)_Real'], label = 'Plune Impingement')
plt.plot(min_['#X'], min_['Displacement(mm)_Real'], label = 'Uniform Pressure (Min Value)')
# plt.plot(max_['#X'], max_['Displacement(mm)_Real'], label = 'Uniform Pressure (Max Value)')
# plt.plot(mean_['#X'], mean_['Displacement(mm)_Real'], label = 'Uniform Pressure (Mean Value)')
plt.legend()
plt.show()