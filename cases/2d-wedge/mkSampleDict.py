import pandas as pd
import os
import sys
import numpy as np

# Create necessary OpenFOAM config files for the required inlet conditions. 
original_stdout = sys.stdout # Save console as an output stream.
r_inlet = 1.005
r_outlet = 20
p = 3

with open("system/sampleDict", "w") as f:
    sys.stdout = f # print to file. 
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
    print()
    print("type sets;")
    print()    
    print("interpolationScheme cellPoint;")
    print()
    print("setFormat csv;")
    print()
    print("sets")
    print("(")

    i = 0
    theta = np.linspace(0, 1.57, 3)
    for curr_theta in theta:
        # print(curr_theta)
        print("     line" + str(i))
        print("     {")
        print("         type uniform;")
        print("         axis distance;")
        x1 = str(round(r_inlet * np.cos(curr_theta), p))
        y1 = str(round(r_inlet * np.sin(curr_theta), p))
        x2 = str(round(r_outlet * np.cos(curr_theta), p))
        y2 = str(round(r_outlet * np.sin(curr_theta), p))
        print("         start ( " + x1 + " " +  y1 + " 0.0 );")
        print("         end ( " + x2 + " " +  y2 + " 0.0 );")
        print("         nPoints 1000;")
        print("     }")
        print()
        i = i + 1
    print(");")
    print()
    print("fields (U_Ar Ttra_Ar rhoN_Ar);")
    print("// ************************************************************************* //")
    f.close()
