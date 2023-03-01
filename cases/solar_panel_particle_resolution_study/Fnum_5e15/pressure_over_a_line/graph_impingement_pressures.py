import pandas as pd
import matplotlib.pyplot as plt

h05 = pd.read_csv('width_05m.csv')
h10 = pd.read_csv('width_10m.csv')
h20 = pd.read_csv('width_20m.csv')


print(h05.dtypes)
print(h05.columns)
# print(plume['#X'])
# plume.plot(kind = 'line', x = '#X', y='Displacement(mm)_Real')
# max_.plot(kind = 'scatter', x = '#X', y='Displacement(mm)_Real')
# mean_.plot(kind = 'scatter', x = '#X', y='Displacement(mm)_Real')
# min_.plot(kind = 'scatter', x = '#X', y='Displacement(mm)_Real')
# # plt.show()
# print(max_['#X'])
# print(mean_['#X'])
# print(min_['#X'])
# print()
# print(type(plume))

# Normalize data
# h05['norm_p_Ar'] = h05['p_Ar'] / max(h05['p_Ar'])
# h10['norm_p_Ar'] = h10['p_Ar'] / max(h05['p_Ar'])
# h20['norm_p_Ar'] = h20['p_Ar'] / max(h05['p_Ar'])

# Convert to mPa
h05['p_Ar'] = 10e9 * h05['p_Ar']
h10['p_Ar'] = 10e9 * h10['p_Ar']
h20['p_Ar'] = 10e9 * h20['p_Ar']


# print(h05['norm_p_Ar'] )

# Label graph with bold characters
font_axis_publish = {
      'color':  'black',
      'weight': 'bold',
      'size': 16,
      }

# # # Graph Results
plt.title('Plume Impingement Pressure (mPa)', fontdict = font_axis_publish)
plt.ylabel('Pressure (mPa)', fontdict = font_axis_publish)
plt.xlabel('Distance Along the Width of Solar Panel (m)', fontdict = font_axis_publish)

plt.plot(h05['Points_1'], h05['p_Ar'], label = 'Plune Impingement (h = 5 m)')
plt.plot(h10['Points_1'], h10['p_Ar'], label = 'Plune Impingement (h = 10 m)')
plt.plot(h20['Points_1'], h20['p_Ar'], label = 'Plune Impingement (h = 10 m)')

# plt.yscale('symlog')
plt.legend()
plt.show()