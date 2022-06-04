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

def readFaceLabels(patch):

    path = 'constant/polyMesh/sets/' + patch
    nFaces = readNfaces(patch)
    face_labels = []

    # print(path)
    with open(path, 'r') as f:
        # print(f.read())

        line = f.readline()
        while str(nFaces) not in line:
            # print(str(nFaces))
            line = f.readline()

        line = f.readline()

        for i in range(nFaces):
            face_labels.append(int(f.readline()))

        f.close()

    return face_labels

def readPointLabels(face_albels):
    path = 'constant/polyMesh/faces'

    # nFaces for the entire mesh.
    tNFaces = 204717
    point_labels = []
    with open(path, 'r') as f:
        line = f.readline()
        while str(tNFaces) not in line:
            # print(str(nFaces))
            line = f.readline()
        print(line)
        line = f.readline()

        for i in range(tNFaces):
            # print(f.readline())
            line = f.readline()
            # Parse line to isolate numbers in the string.
            line = line.split('(')
            line = line[1].split(' ')
            line[-1] = line[-1].split(')')[0]
            print(line)
            point_labels.append(line)

    return point_labels
        # print(line)
        # print('printing line',line)/
# Reading labels of faces.
# # with open('constant/polyMesh/sets/inflow', 'r') as f:
# #     print(f.read())
patch = 'inflow'
print(readNfaces(patch))
print(type(readNfaces('inflow')))

face_labels = readFaceLabels(patch)
# print(face_labels)
readPointLabels(face_labels)