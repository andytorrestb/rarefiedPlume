import os
import sys

def get_parent_cases(curr_dir_path):
    cases = os.listdir(curr_dir_path + '/Cases')
    pop = cases.index('baseCase')
    cases.pop(pop)
    return cases

def mkSampleDict(sampleDict_path):
    original_stdout = sys.stdout # Save console as an output stream.

    with open(sampleDict_path, "w") as f:
        print("/*--------------------------------*- C++ -*----------------------------------*\\")
        print("| =========                 |                                                 |")
        print("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |")
        print("|  \\\\    /   O peration     | Version:  v1706                                 |")
        print("|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |")
        print("|    \\\\/     M anipulation  |                                                 |")
        print("\*---------------------------------------------------------------------------*/")
        print("FoamFile")
        print("{")
        print("     version     2.0;")
        print("     format      ascii;")
        print("     class       dictionary;")
        print("     location    \"system\";")
        print("     object      sampleDict;")
        print("}")
        print("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")

        print("type sets;")

        print("interpolationScheme cellPoint;")
        print("setFormat csv;")
        print("sets")
        print("(")
        print("     radialLine")
        print("     {")
        print("             type uniform;")
        print("             axis distance;")
        print("             start (0.74246 0.74246 0.0;")
        print("             end   (14.1421 14.1421 0.0);")
        print("             nPoints 100;")
        print("     }")
        print(");")

        print()

        print("fields (U_Ar Ttra_ar rhoN_ar);")
        print("// ************************************************************************* //")
        sys.stdout = original_stdout
        f.close()

    return

def main():
    # Get working directory
    curr_dir_path = os.path.dirname(os.path.realpath(__file__))
    parent_cases = get_parent_cases(curr_dir_path)

    # Loop through every produced case and sample data along the radial line. 
    for case_parent in parent_cases:
        for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
            curr_case = curr_dir_path + '/Cases/' + case_parent + '/' + case_child
            print(curr_case)

            sampleDict_path =  curr_case + "/system/sampleDict"
            mkSampleDict(sampleDict_path)

            os.system("cat " + curr_case + "/system/sampleDict")
            input()

main()