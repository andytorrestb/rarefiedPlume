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
        # os.system("touch " + curr_case + "/0/boundaryT")
        # os.system("touch " + curr_case + "/0/boundaryU")
        # os.system("touch " + curr_case + "/0/boundaryNumberDensity_Ar")
        os.system("python3 " + curr_case + '/processInflowData.py 743.397 802.2867 6.0201e20')
        os.system("decomposePar -case " + curr_case)
        os.system("mpirun -np 8 dsmcFoam+ -case " + curr_case)
        # os.system("dsmcFoam+ -case " + curr_case)
   # print(case_parent + " : " + str(os.listdir(dir_path + '/Cases/' + case)))