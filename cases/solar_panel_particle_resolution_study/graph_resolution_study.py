import pandas as pd
import matplotlib.pyplot as plt

# Load sampled data into a dictionary. The key is the parameter changed. 
cases = {}
file_name = 'line_h_7-5_offset_1-1.csv'
cases['2.5e15'] = pd.read_csv('Fnum_2.5e15/' + file_name)
cases['5e15'] = pd.read_csv('Fnum_5e15/' + file_name)
cases['1e16'] = pd.read_csv('Fnum_1e16/' + file_name)

# Plot data.
for case in cases:
    plt.plot(cases[case]['Points_1'], cases[case]['p_Ar'], label = 'Fnum = ' + str(case))    

# Label graph with bold characters.
font_axis_publish = {
      'color':  'black',
      'weight': 'bold',
      'size': 22,
      }

# Set Title and Label Text.
plt.title('Plume Impingement Pressure (mPa)', fontdict = font_axis_publish)
plt.xlabel('Distance Along Width of the Solar Panel (m)', fontdict = font_axis_publish)
plt.ylabel('Plume Impingement Pressure (mPa)', fontdict = font_axis_publish)

# Do the thing. 
plt.legend()
plt.show()

