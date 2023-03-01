import pandas as pd
import matplotlib.pyplot as plt

h05 = pd.read_csv('line_h_5_offset_1-1.csv')
h075 = pd.read_csv('line_h_7-5_offset_1-1.csv')
h10 = pd.read_csv('line_h_10_offset_1-1.csv')
h20 = pd.read_csv('line_h_20_offset_1-1.csv')

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
factor = 1
# h05['p_Ar'] = factor * h05['p_Ar']
# h075['p_Ar'] = factor * h075['p_Ar']
# h10['p_Ar'] = factor * h10['p_Ar']
h20['p_Ar'] = factor * h20['p_Ar']


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

factor = 1000
# plt.plot(h05['Points_1'], factor*h05['p_Ar'], label = 'Plune Impingement (h = 5 m)')
# plt.plot(h075['Points_1'], factor*h075['p_Ar'], label = 'Plune Impingement (h = 7.5 m)')
# plt.plot(h10['Points_1'], factor*h10['p_Ar'], label = 'Plune Impingement (h = 10 m)')
plt.plot(h20['Points_1'], factor*h20['p_Ar'], label = 'Plune Impingement (h = 20 m)', color = 'red')

# plt.yscale('symlog')
plt.legend()
plt.show()