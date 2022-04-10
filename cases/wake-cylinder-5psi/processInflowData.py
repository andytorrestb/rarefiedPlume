# Author: Andy Torres
# Purpose: Load inflow data from .dat file and write it to the proper boundary dictoinary in the 0/ (initial conditions) folder.

import pandas as pd
import os
import sys

# Load data, assign columns
df = pd.read_table("inflowData.dat", sep = "\s+")
df.columns = ["x coordinate", "x velocity", "y velocity", "temperature", "Nd"]

# Add y - coordinates to data set.

# Equation of a circle
# (x-h)^2 + (y-k)^2 = r^2

# Solving for y = f(x) yields,
# y = k + [r^2 -(x-h)^2] ^ 0.5

k = 0
h = -0.37465 # 11.75 + 3 inches converted to meters, in negative direction
r = 0.1524 # 6 inches converted to meters
df["y coordinate"] = k + pow(r**2 -pow(df["x coordinate"]-h, 2), 0.5)

# Re- order columns
df = df[["x coordinate", "y coordinate", "x velocity", "y velocity", "temperature", "Nd"]]
print(df.shape)

velo_data = df[["x velocity", "y velocity"]]

print(velo_data.shape[0])


original_stdout = sys.stdout

# =============== Printing boundaryU =======================================
with open("0/boundaryU", 'w') as f:
    sys.stdout = f
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
    print("     class       volVectorField;")
    print("     location    \"0\";")
    print("     object      boundaryU;")
    print("}")
    print("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")

    print("dimensions      [0 1 -1 0 0 0 0];")


    print("internalField   uniform (0 0 0);")

    print("boundaryField")
    print("{")
    print("    inflow")
    print("    {")
    print("        type            fixedValue;")
    print("        value           nonuniform List<vector>") 
    print("        " + str(velo_data.shape[0]))
    print("        (")

    for row in velo_data.iterrows():
        # print(type(row))
        # print(type(row[0]), type(row[1]))
        data = row[1]
        # print(data.name, data[0], data[1])
        print("            (" + str(data[0]) + " " + str(data[1]) + " 0.0)")
        # print()
    print("        );")
    print("    }")
    print("    vacuum")
    print("    {")
    print("         type zeroGradient;")
    print("    }")
    print()
    print("    noSolution")
    print("    {")
    print("         type empty;")
    print("    }")
    print()
    print("    sym")
    print("    {")
    print("         type symmetry;")
    print("    }")
    print()
    print("    cylinder")
    print("    {")
    print("         type  calculated;")
    print("         value uniform (0 0 0);")
    print("    }")
    print()
    print("    plate")
    print("    {")
    print("         type  calculated;")
    print("         value uniform (0 0 0);")
    print("    }")
    print("}")
    sys.stdout = original_stdout
    f.close()

# =============== Printing boundaryT =======================================
with open("0/boundaryT", 'w') as f:
    sys.stdout = f
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
    print("     class       volVectorField;")
    print("     location    \"0\";")
    print("     object      boundaryT;")
    print("}")
    print("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")

    print("dimensions      [0 0 0 1 0 0 0];")


    print("internalField   uniform (0 0 0);")

    print("boundaryField")
    print("{")
    print("    inflow")
    print("    {")
    print("        type            fixedValue;")
    print("        value           nonuniform List<vector>") 
    print("        " + str(velo_data.shape[0]))
    print("        (")

    for row in df.iterrows():
        # print(type(row))
        # print(row)
        # print(row[1][4])
        # print(type(row[0]), type(row[1]))
        # data = row[1]
        # print(data.name, data[0], data[1])
        print("            ( " + str(row[1][4]) + " 0.0 0.0 )")
        # print()
    print("        );")
    print("    }")
    print("    vacuum")
    print("    {")
    print("         type zeroGradient;")
    print("    }")
    print()
    print("    noSolution")
    print("    {")
    print("         type empty;")
    print("    }")
    print()
    print("    sym")
    print("    {")
    print("         type symmetry;")
    print("    }")
    print()
    print("    cylinder")
    print("    {")
    print("         type  calculated;")
    print("         value uniform (0 0 0);")
    print("    }")
    print()
    print("    plate")
    print("    {")
    print("         type  calculated;")
    print("         value uniform (0 0 0);")
    print("    }")
    print("}")
    sys.stdout = original_stdout
    f.close()

# =============== Printing boundaryNumberDensity_Ar =======================================
with open("0/boundaryNumberDensity_Ar", 'w') as f:
    sys.stdout = f
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
    print("     class       volScalarField;")
    print("     location    \"0\";")
    print("     object      boundaryNumberDensity_Ar;")
    print("}")
    print("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")

    print("dimensions      [0 -3 0 0 0 0 0];")


    print("internalField   uniform 0;")

    print("boundaryField")
    print("{")
    print("    inflow")
    print("    {")
    print("        type            fixedValue;")
    print("        value           nonuniform List<scalar>") 
    print("        " + str(velo_data.shape[0]))
    print("        (")

    for row in df.iterrows():
        # print(type(row))
        # print(row)
        # print(row[1][4])
        # print(type(row[0]), type(row[1]))
        # data = row[1]
        # print(data.name, data[0], data[1])
        print("              " + str(row[1][5]))
        # print()
    print("        );")
    print("    }")
    print("    vacuum")
    print("    {")
    print("         type zeroGradient;")
    print("    }")
    print()
    print("    noSolution")
    print("    {")
    print("         type empty;")
    print("    }")
    print()
    print("    sym")
    print("    {")
    print("         type symmetry;")
    print("    }")
    print()
    print("    cylinder")
    print("    {")
    print("         type  calculated;")
    print("         value uniform 0;")
    print("    }")
    print()
    print("    plate")
    print("    {")
    print("         type  calculated;")
    print("         value uniform 0;")
    print("    }")
    print("}")
    sys.stdout = original_stdout
    f.close()