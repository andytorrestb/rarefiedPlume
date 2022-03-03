import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# Find path for cases
curr_dir_path = os.path.dirname(os.path.realpath(__file__))
cases = os.listdir(curr_dir_path + '/Cases')
pop = cases.index('baseCase')
cases.pop(pop)

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


# Traverse through case structure, and graph all relevant data. 
for case_parent in cases:
    parent_path = curr_dir_path + '/Cases/' + case_parent
    for case_child in os.listdir(parent_path):

        curr_case = parent_path + '/' + case_child
        print(curr_case)
        sim = pd.read_csv(
            curr_case + '/results.csv'
        )

        print(sim['Ttra_Ar'])

        # print(sim.columns)

        sim = sim[['Points_Magnitude', 'rhoN_Ar', 'Ttra_Ar']].dropna()
        sim['rhoN_Ar'] = sim['rhoN_Ar'] / 8.773e20
        sim['Ttra_Ar'] = sim['Ttra_Ar'] / 1e3
         

        plt.title('Analytical vs Digitized Data', fontdict = font_axis_publish)
        plt.ylabel('n/n*', fontdict = font_axis_publish)
        plt.xlabel('r', fontdict = font_axis_publish)

        plt.plot(digi_n['r'], digi_n['n_nstar'], label = 'digitized')
        plt.plot(sim['Points_Magnitude'], sim['rhoN_Ar'], label = 'dsmc')
        plt.legend()
        plt.savefig(curr_case + '/digitized_vs_analytical_n.png')
        plt.close()

        plt.ylabel('T/T*', fontdict = font_axis_publish)
        plt.xlabel('r', fontdict = font_axis_publish)

        plt.plot(digi_T['r'], digi_T['T_Tstar'], label = 'digitized')
        plt.plot(sim['Points_Magnitude'], sim['Ttra_Ar'], label = 'dsmc')
        plt.legend()
        plt.savefig(curr_case + '/digitized_vs_analytical_T.png')
        plt.close()
