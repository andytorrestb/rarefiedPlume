import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def processSampledData(df):
    df = df[['distance', 'rhoN_Ar', 'Ttra_Ar']].dropna()
    df['distance'] = df['distance'] + 1.005
    df['rhoN_Ar'] = df['rhoN_Ar'] / 8.377e20
    df['Ttra_Ar'] = df['Ttra_Ar'] / 1000.0
    return df

# Find path for cases
curr_dir_path = os.path.dirname(os.path.realpath(__file__))

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
    curr_dir_path + '/T_Tstar_radius_DAC.dat',
    header = 0,
    sep = '\t',
    names = ['r', 'T_Tstar']
    )

# Read in simulated data. 
sim = []
for i in range(3):
    path_to_data = '/postProcessing/sampleDict/0.1/line' + str(i) + '_Ttra_Ar_rhoN_Ar.csv'
    print(path_to_data)
    sim.append(processSampledData(pd.read_csv(curr_dir_path + path_to_data)))


# Produce analytical data
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

Ma_domain = np.linspace(1, 25, 100) 
ga = 1.67
TTt = TTt_Ma(Ma_domain, ga = ga)
rrt = rrt_Ma(Ma_domain, ga = ga)
nnt = nnt_Ma(Ma_domain, ga = ga)


# Report # and angle for each line
i = 0
theta = 90.0 / 3.0
degree_sign = u'\N{DEGREE SIGN}'

# Plot simulated data
plt.title('OpenFOAM Samples vs DAC', fontdict = font_axis_publish)
plt.ylabel('T/T**', fontdict = font_axis_publish)
plt.xlabel('Radial distance, r (m)', fontdict = font_axis_publish)
plt.plot(digi_T['r'], digi_T['T_Tstar'], label = 'DAC')
plt.plot(rrt, TTt, label = 'Analytical')
for line in sim:
    print(line)
    label = 'OpenFOAM line' + str(i) + ' theta ' + str(i*theta) + degree_sign
    i = i + 1
    plt.plot(line['distance'], line['Ttra_Ar'], label = label)
    # plt.plot(line['distance'], line['Ttra_Ar'])
    

plt.legend()
plt.yscale('log')
plt.ylim(bottom = 1e-4, top = 1)
plt.savefig(curr_dir_path + '/DAC_vs_OF_samples_T.png')
plt.close()