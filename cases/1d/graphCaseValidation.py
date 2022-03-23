import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# Find path for cases
curr_dir_path = os.path.dirname(os.path.realpath(__file__))
# print(curr_dir_path)
# cases = os.listdir(curr_dir_path + '/Cases')
# pop = cases.index('baseCase')
# cases.pop(pop)

# Label graph with bold characters
font_axis_publish = {
      'color':  'black',
      'weight': 'bold',
      'size': 22,
      }

# Read in digitized data
digi_n = pd.read_csv(
    curr_dir_path + '/n_nstar_radius.dat',
    header = 0,
    sep = '\t',
    names = ['r', 'n_nstar']
    )

digi_T = pd.read_csv(
    curr_dir_path + '/T_Tstar_radius.dat',
    header = 0,
    sep = '\t',
    names = ['r', 'T_Tstar']
    )

sim = pd.read_csv(
    curr_dir_path + '/postProcessing/sampleDict/0.36/horizontalLine_Ttra_Ar_rhoN_Ar.csv'
    )

print(sim['Ttra_Ar'])

sim = sim[['x', 'rhoN_Ar', 'Ttra_Ar']].dropna()
sim['rhoN_Ar'] = sim['rhoN_Ar'] / 6.02e20
sim['Ttra_Ar'] = sim['Ttra_Ar'] / 800.0
    

plt.title('DSMC vs DAC', fontdict = font_axis_publish)
plt.ylabel('n/n*', fontdict = font_axis_publish)
plt.xlabel('r', fontdict = font_axis_publish)

plt.plot(digi_n['r'], digi_n['n_nstar'], label = 'digitized (DAC)')
plt.plot(sim['x'], sim['rhoN_Ar'], label = 'simulated (DSMC)')
plt.legend()
plt.savefig(curr_dir_path + '/digitized_vs_analytical_n.png')
plt.close()

plt.title('DSMC vs DAC', fontdict = font_axis_publish)
plt.ylabel('T', fontdict = font_axis_publish)
plt.xlabel('r', fontdict = font_axis_publish)

plt.plot(digi_T['r'], digi_T['T_Tstar'], label = 'digitized (DAC)')
plt.plot(sim['x'], sim['Ttra_Ar'], label = 'simulated (DSMC)')
plt.legend()
plt.savefig(curr_dir_path + '/digitized_vs_analytical_T.png')
plt.close()
