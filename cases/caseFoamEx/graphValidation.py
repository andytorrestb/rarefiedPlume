import os
import pandas as pd
import matplotlib.pyplot as plt

# Get working directory
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

def oneGraph():
    # Loop through every produced case and sample data along the radial line.
    for case_parent in cases:
        for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
            curr_case = curr_dir_path + '/Cases/' + case_parent + '/' + case_child
            os.system("postProcess -case " + curr_case + " -func sampleDict")

    # Loop through cases, plotting n/n* on the same graph
    plt.title('DSMC vs DAC', fontdict = font_axis_publish)
    plt.ylabel('n/n*', fontdict = font_axis_publish)
    plt.xlabel('r', fontdict = font_axis_publish)
    plt.plot(digi_n['r'], digi_n['n_nstar'], label = 'digitized (DAC)')
    i = 0

    for case_parent in cases:
        for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
            sim = pd.read_csv(
                curr_case + '/postProcessing/sampleDict/0.15/radialLine_Ttra_Ar_rhoN_Ar.csv'
            )
            print(sim)
            sim = sim[['distance', 'rhoN_Ar', 'Ttra_Ar']].dropna()
            sim['rhoN_Ar'] = sim['rhoN_Ar'] / 6.02e20
            sim['Ttra_Ar'] = sim['Ttra_Ar'] / 800.0
            plt.plot(sim['distance'], sim['rhoN_Ar'], label = 'run #' + str(i) + ' (DSMC)')
            i = i + 1

    plt.legend()
    plt.yscale('log')
    plt.savefig(curr_dir_path + '/digitized_vs_analytical_n.png')
    plt.close()

    # Loop through cases, plotting T/T* on the same graph
    plt.title('DSMC vs DAC', fontdict = font_axis_publish)
    plt.ylabel('T/T*', fontdict = font_axis_publish)
    plt.xlabel('r', fontdict = font_axis_publish)
    plt.plot(digi_T['r'], digi_T['T_Tstar'], label = 'digitized (DAC)')
    i = 0
    for case_parent in cases:
        for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
            sim = pd.read_csv(
                curr_case + '/postProcessing/sampleDict/0.15/radialLine_Ttra_Ar_rhoN_Ar.csv'
            )
            sim = sim[['distance', 'rhoN_Ar', 'Ttra_Ar']].dropna()
            sim['rhoN_Ar'] = sim['rhoN_Ar'] / 6.02e20
            sim['Ttra_Ar'] = sim['Ttra_Ar'] / 800.0
            plt.plot(sim['distance'], sim['Ttra_Ar'], label = 'run #' + str(i) + ' (DSMC)')
            i = i + 1

    plt.legend()
    plt.yscale('log')
    plt.savefig(curr_dir_path + '/digitized_vs_analytical_T.png')
    plt.close()

def nGraphs():
    # Loop through every produced case and sample data along the radial line.
    for case_parent in cases:
        print(case_parent)
        for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
            curr_case = curr_dir_path + '/Cases/' + case_parent + '/' + case_child
            print(curr_case)

oneGraph()
# nGraphs()