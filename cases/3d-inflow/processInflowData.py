def readNfaces(patch):

    path = 'constant/polyMesh/boundary'

    with open(path, 'r') as f:
        line = f.readline()

        # finds the given string 'patch', 
        while patch not in line:
            line = f.readline()
        # skips a few lines 
        line = f.readline()
        line = f.readline()
        line = f.readline()
        f.close()
        # parses the string to return nFaces as an int.
        return int(line.split(' ', 14)[-1][0:-2])

    return

# Reading labels of faces.
# # with open('constant/polyMesh/sets/inflow', 'r') as f:
# #     print(f.read())

print(readNfaces('inflow'))
print(type(readNfaces('inflow')))