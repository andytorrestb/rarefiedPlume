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

# Read in simulated data. 
sim = pd.read_csv(
    curr_dir_path + '/postProcessing/sampleDict/0.3/horizontalLine_Ttra_Ar_rhoN_Ar.csv'
    )

# Used to see what the values trend to. 
print(sim['Ttra_Ar'])

sim = sim[['x', 'rhoN_Ar', 'Ttra_Ar']].dropna()
sim['rhoN_Ar'] = sim['rhoN_Ar'] / 6.02e20
sim['Ttra_Ar'] = sim['Ttra_Ar'] / 800.0
    

# Producde Analytical Data
def TTt_Ma(Ma, ga = 1.4):
    return (ga + 1) / (2 + (ga - 1) * Ma ** 2)

def rrt_Ma(Ma, ga = 1.4):
  rrt = (1 / TTt_Ma(Ma, ga)) ** ((ga + 1) / (ga - 1))
  rrt = np.sqrt(np.sqrt(rrt) / Ma)
  return rrt

def nnt_Ma(Ma, ga = 1.4):
  return TTt_Ma(Ma, ga) ** (1 / (ga - 1))

def a(T, ga = 1.4, R = 287):
    return np.sqrt(ga * R * T)

Ma_domain = np.linspace(1, 10.5, 100)
TTt = TTt_Ma(Ma_domain)
rrt = rrt_Ma(Ma_domain)
nnt = nnt_Ma(Ma_domain)

print("Printing rrt")
print(rrt)

# Graph Results
plt.title('DSMC vs DAC', fontdict = font_axis_publish)
plt.ylabel('n/n*', fontdict = font_axis_publish)
plt.xlabel('r', fontdict = font_axis_publish)

plt.plot(digi_n['r'], digi_n['n_nstar'], label = 'digitized (DAC)')
plt.plot(sim['x'], sim['rhoN_Ar'], label = 'simulated (DSMC)')
plt.plot(rrt, nnt, label = 'analytical (Isentropic Expansion)')
plt.legend()
plt.yscale('log')
plt.savefig(curr_dir_path + '/digitized_vs_analytical_n.png')
plt.close()

plt.title('DSMC vs DAC', fontdict = font_axis_publish)
plt.ylabel('T/T*', fontdict = font_axis_publish)
plt.xlabel('r', fontdict = font_axis_publish)

plt.plot(digi_T['r'], digi_T['T_Tstar'], label = 'digitized (DAC)')
plt.plot(sim['x'], sim['Ttra_Ar'], label = 'simulated (DSMC)')
plt.plot(rrt, TTt, label = 'analytical (Isentropic Expansion)')
plt.legend()
plt.yscale('log')
plt.savefig(curr_dir_path + '/digitized_vs_analytical_T.png')
plt.close()
