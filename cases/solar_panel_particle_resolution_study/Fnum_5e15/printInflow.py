import sys
import os

def printInflowSurface(inflow):
    face_rhoN = inflow[0]
    face_U = inflow[1]
    face_T = inflow[2]


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
        print("        " + str(len(face_U)))
        print("        (")

        for face in face_U:
            # print(type(row))
            # print(type(row[0]), type(row[1]))
            data = face_U[face]
            # print(data[3], data[4])
            print("            (" + str(data[0]) + " " + str(data[1]) + " " + str(data[2]) + ")")
            # print()
        print("        );")
        print("    }")
        print("    vacuum")
        print("    {")
        print("         type zeroGradient;")
        print("    }")
        print()
        print("    sym")
        print("    {")
        print("         type symmetry;")
        print("    }")
        # print()
        # print("    cylinder")
        # print("    {")
        # print("         type  calculated;")
        # print("         value uniform (0 0 0);")
        # print("    }")
        print()
        print("    panel")
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
        print("        " + str(len(face_T)))
        print("        (")

        for face in face_T:
            # print(type(row))
            # print(row)
            # print(row[1][4])
            # print(type(row[0]), type(row[1]))
            data = face_T[face]
            # print(data.name, data[0], data[1])
            print("            ( " + str(data) + " 0.0 0.0 )")
            # print("            " + str(data[5]))
            # print()
        print("        );")
        print("    }")
        print("    vacuum")
        print("    {")
        print("         type zeroGradient;")
        print("    }")
        print()
        # print("    sym")
        # print("    {")
        # print("         type symmetry;")
        # print("    }")
        # print()
        # print("    cylinder")
        # print("    {")
        # print("         type  calculated;")
        # print("         value uniform (0 0 0);")
        # print("    }")
        # print()
        # print("    panel")
        # print("    {")
        # print("         type  calculated;")
        # print("         value uniform (0 0 0);")
        # print("    }")
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
        print("        " + str(len(face_rhoN)))
        print("        (")

        for face in face_rhoN:
            # print(type(row))
            # print(row)
            # print(row[1][4])
            # print(type(row[0]), type(row[1]))
            data = face_rhoN[face]
            # print(data.name, data[0], data[1])
            print("              " + str(data))
            # print()
        print("        );")
        print("    }")
        print("    vacuum")
        print("    {")
        print("         type zeroGradient;")
        print("    }")
        print()
        print("    sym")
        print("    {")
        print("         type symmetry;")
        print("    }")
        print()
        # print("    cylinder")
        # print("    {")
        # print("         type  calculated;")
        # print("         value uniform 0;")
        # print("    }")
        # print()
        print("    panel")
        print("    {")
        print("         type  calculated;")
        print("         value uniform 0;")
        print("    }")
        print("}")
        sys.stdout = original_stdout
        f.close()