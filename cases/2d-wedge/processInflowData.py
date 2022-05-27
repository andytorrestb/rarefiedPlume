# Author: Andy Torres
# Purpose: Produce inflow data of constant value in the radial direction according to a given magintude. 

import pandas as pd
import os
import sys
import numpy as np

# Declare input data, angles from 0 to 90 degrees
theta = np.linspace(0, 1.57, 33)
# print(theta)

p = 4 # Precision for rounding
radius = 1.05 # inflow radius

# Check for proper arguments before processing
if (len(sys.argv) < 3):
    print("Please enter 3 arguments")
    print("processInflowData.py U_mag T rhoN")
    sys.exit()

# Co-ordinates according to range of theta
x = np.round(radius * np.cos(theta), p)
print(x)
y = np.round(radius * np.sin(theta), p)
print(y)

# Velocity components accoring to user input and range of theta
velocity = float(sys.argv[1])
u = np.round(velocity* np.cos(theta), p)
print(u)
v = np.round(velocity* np.sin(theta), p)
print(v)

# Temperature Values according to user input. 
temperature = float(sys.argv[2])
Tx = np.round(temperature * np.cos(0), p)
Ty = np.round(temperature * np.sin(0), p)
print(Tx)
print(Ty)

rhoN = float(sys.argv[3])
Nd = np.linspace(rhoN, rhoN, 33)
print(Nd)

# Save produced data in a dataframe (table of values)
velo_data = pd.DataFrame(columns = ['theta','x', 'y', 'u', 'v', 'Tx', 'Ty', 'Nd'])
velo_data['theta'] = theta
velo_data['x'] = x
velo_data['y'] = y
velo_data['u'] = u
velo_data['v'] = v
velo_data['Tx'] = Tx
velo_data['Ty'] = Ty
velo_data['Nd'] = Nd

# Print to confirm values are reasonable.
print(velo_data)

# Create necessary OpenFOAM config files for the required inlet conditions. 
original_stdout = sys.stdout # Save console as an output stream.

# =============== Printing boundaryU =======================================
with open("0/boundaryU", 'w') as f:
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
    print("        " + str(u.shape[0]))
    print("        (")

    for row in velo_data.iterrows():
        # print(type(row))
        # print(type(row[0]), type(row[1]))
        data = row[1]
        # print(data[3], data[4])
        print("            (" + str(data[3]) + " " + str(data[4]) + " 0.0)")
        # print()
    print("        );")
    print("    }")
    print("    vacuum")
    print("    {")
    print("         type zeroGradient;")
    print("    }")
    print()
    print("    front")
    print("    {")
    print("         type symmetry;")
    print("    }")
    print()
    print("    back")
    print("    {")
    print("         type symmetry;")
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

    for row in velo_data.iterrows():
        # print(type(row))
        # print(row)
        # print(row[1][4])
        # print(type(row[0]), type(row[1]))
        data = row[1]
        # print(data.name, data[0], data[1])
        print("            ( " + str(data[5]) + " " + str(data[6]) + " 0.0 )")
        # print()
    print("        );")
    print("    }")
    print("    vacuum")
    print("    {")
    print("         type zeroGradient;")
    print("    }")
    print()
    print("    front")
    print("    {")
    print("         type symmetry;")
    print("    }")
    print()
    print("    back")
    print("    {")
    print("         type symmetry;")
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

    for row in velo_data.iterrows():
        # print(type(row))
        # print(row)
        # print(row[1][4])
        # print(type(row[0]), type(row[1]))
        data = row[1]
        # print(data.name, data[0], data[1])
        print("              " + str(data[7]))
        # print()
    print("        );")
    print("    }")
    print("    vacuum")
    print("    {")
    print("         type zeroGradient;")
    print("    }")
    print()
    print("    front")
    print("    {")
    print("         type symmetry;")
    print("    }")
    print()
    print("    back")
    print("    {")
    print("         type symmetry;")
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