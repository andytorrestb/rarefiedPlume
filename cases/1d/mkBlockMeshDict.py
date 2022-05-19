import os
import sys
import math

def printHeader():
    print("/*--------------------------------*- C++ -*----------------------------------*\\")
    print("| =========                 |                                                 |")
    print("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |")
    print("|  \\\\    /   O peration     | Version:  v1706                                 |")
    print("|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |")
    print("|    \\\\/     M anipulation  |                                                 |")
    print("\*---------------------------------------------------------------------------*/")

def printConfig():
    print("FoamFile")
    print("{")
    print("     version     2.0;")
    print("     format      ascii;")
    print("     class       dictionary;")
    print("     object      blockMeshDict;")
    print("}")
    print("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")   
    print()

def printScale():
    print("scale        1;")
    print()

def printPoints():
    r_inlet = float(sys.argv[1])
    r_outlet = float(sys.argv[2])
    alpha = (math.pi / 180) * float(sys.argv[3]) 

    p = 3 # precision

    x0 = str(0)
    x1 = str(round(r_inlet * math.cos(alpha), p))
    x2 = str(round(r_outlet * math.cos(alpha), p))

    y0 = str(round(-1 * r_outlet * math.sin(alpha), p))
    y1 = str(round(-1 * r_inlet * math.sin(alpha), p))
    y2 = str(round(r_inlet * math.sin(alpha), p))
    y3 = str(round(r_outlet * math.sin(alpha), p))

    z0 = y0
    z1 = y1
    z2 = y2
    z3 = y3

    p0 = x1 + " " + y1 + " " + z1
    p1 = x2 + " " + y0 + " " + z0
    p2 = x2 + " " + y3 + " " + z0
    p3 = x1 + " " + y2 + " " + z1

    p4 = x1 + " " + y1 + " " + z2
    p5 = x2 + " " + y0 + " " + z3
    p6 = x2 + " " + y3 + " " + z3
    p7 = x1 + " " + y2 + " " + z2

    print("vertices")
    print("(")
    print("    (" + p0 + ") // 0")
    print("    (" + p1 + ") // 1")
    print("    (" + p2 + ") // 2")
    print("    (" + p3 + ") // 3")
    print("    (" + p4 + ") // 4")
    print("    (" + p5 + ") // 5")
    print("    (" + p6 + ") // 6")
    print("    (" + p7 + ") // 7")
    print(");")
    print()

def printBlock():
    print("blocks")
    print("(")
    print("     hex (0 1 2 3 4 5 6 7) (255 1 1) simpleGrading (1 1 1)")
    print(");")
    print()

def printEdges():
    print("edges")
    print("(")
    print(");")
    print()

def printBoundary():
    print("patches")
    print("(")
    print(");")
    print()

def printMergePatchPairs():
    print("mergePatchPairs")
    print("(")
    print(");")

    
# Find path for cases
curr_dir_path = os.path.dirname(os.path.realpath(__file__))

print(curr_dir_path)

if (len(sys.argv) < 3):
    print("Please enter 3 arguments")
    print("mkBlockMeshDict.py r_inlet r_outlet alpha")
    sys.exit()

# Create necessary OpenFOAM config file for running blockMesh
original_stdout = sys.stdout # Save console as an output stream.

with open("system/blockMeshDict", 'w') as f:
    sys.stdout = f
    printHeader()
    printConfig()
    printScale()
    printPoints()
    printBlock()
    printEdges()
    printBoundary()
    print("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")

    sys.stdout = original_stdout
    f.close()