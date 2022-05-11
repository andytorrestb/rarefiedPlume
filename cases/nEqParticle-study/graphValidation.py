import os

curr_dir_path = os.path.dirname(os.path.realpath(__file__))
cases = os.listdir(curr_dir_path + '/Cases')
pop = cases.index('baseCase')
cases.pop(pop)

# Loop through every produced case. 
for case_parent in cases:
    for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
        curr_case = curr_dir_path + '/Cases/' + case_parent + '/' + case_child
        os.system("postProcess -case " + curr_case + " -func sampleDict")
        print(curr_case)