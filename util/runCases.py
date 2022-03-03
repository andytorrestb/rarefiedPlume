import os

# Find path for cases
curr_dir_path = os.path.dirname(os.path.realpath(__file__))
cases = os.listdir(curr_dir_path + '/Cases')
pop = cases.index('baseCase')
cases.pop(pop)

# Initialise and run solver
for case_parent in cases:
    for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
        curr_case = curr_dir_path + '/Cases/' + case_parent + '/' + case_child
        os.system("dsmcInitialise+ -case " + curr_case)
        os.system("dsmcFoam+ -case " + curr_case)
   # print(case_parent + " : " + str(os.listdir(dir_path + '/Cases/' + case)))