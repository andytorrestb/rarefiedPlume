from meshStats import readMeshStats

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
    # Read labels to all faces in the mesh
    path = 'constant/polyMesh/sets/' + patch
    nFaces = readNfaces(patch)
    face_labels = {}

    # print(path)
    with open(path, 'r') as f:
        # print(f.read())

        line = f.readline()
        while str(nFaces) not in line:
            # print(str(nFaces))
            line = f.readline()

        line = f.readline()

        for i in range(nFaces):
            # face_labels.append(int(f.readline()))
            face_labels[str(i)] = int(f.readline())
        f.close()

    return face_labels

def readFacePoints(face_labels, mesh_stats):
# Read labels for all points in the mesh.
    path = 'constant/polyMesh/faces'
    face_dict = {}
    face_points = {}
    with open(path, 'r') as f:
        line = f.readline()
        while str(mesh_stats['faces']) not in line:
            line = f.readline()
        f.readline()
        # Save all points in mesh to a dictionary.
        for i in range(mesh_stats['faces']):
            # print(i, f.readline())
            face_dict[str(i)] = f.readline()


        # print(len(face_dict))
        # print(len(face_labels))

        for face in face_labels:
            # print(face)
            # print(face, face_dict[str(face)])
            # print(face_dict[str(face)].split("(")[1].split(")")[0].split(' '))
            face_points[str(face)] = face_dict[str(face)].split("(")[1].split(")")[0].split(' ')

    return face_points
    # for i in range(len(face_labels)):
    #     print(i, face_labels[i])


def readPointLabels(tNFaces):
    path = 'constant/polyMesh/faces'

    # nFaces for the entire mesh.
    # tNFaces = 204717
    point_labels = {}
    with open(path, 'r') as f:
        line = f.readline()
        while str(tNFaces) not in line:
            # print(str(nFaces))
            line = f.readline()
        # print(line)
        line = f.readline()

        for i in range(tNFaces):
            # print(f.readline())
            line = f.readline()
            # Parse line to isolate numbers in the string.
            line = line.split('(')
            line = line[1].split(' ')
            line[-1] = line[-1].split(')')[0]
            # print(line)
            point_labels[str(i)] = line

    return point_labels

def processPointCoordinates(nPoints):
    path = 'constant/polyMesh/points'
    point_cooordinates = {}

    with open(path, 'r') as f:
        line = f.readline()
        while str(nPoints) not in line:
            line = f.readline()

        line = f.readline()

        for i in range(nPoints):
            line = f.readline().strip("(").strip('\n').strip(")").split(" ")
            # print(line)
            point_cooordinates[str(i)] = line

    return point_cooordinates

def processCentroids(face_labels, point_labels, point_cooordinates):

    centroid = {}
    for face in face_labels:

        face = str(face)
        print('face', face)
        print('face_points[face]', face_points[face])
        print('point_labels[face]', point_labels[face])
        # p

        print(len(face_points[face]))

        # for i in range(len(face_points[face])):
        #     # print(point_cooordinates[str(face)])
        #     coordinate = str(face_points[face][i])
        #     print(coordinate, point_cooordinates[coordinate])

        print()

        centroid[face] = calculateCentroid(face, face_points, point_cooordinates)
    return centroid

def calculateCentroid(face, face_points, point_cooordinates):

    x_c = 0
    y_c = 0
    z_c = 0

    for i in range(len(face_points[face])):
        # print(point_cooordinates[str(face)])
        coordinate = str(face_points[face][i])
        print(coordinate, point_cooordinates[coordinate])
        x_c = x_c + float(point_cooordinates[coordinate][0])
        y_c = y_c + float(point_cooordinates[coordinate][1])
        z_c = z_c + float(point_cooordinates[coordinate][2])


    x_c = x_c / 3.0
    y_c = y_c / 3.0
    z_c = z_c / 3.0

    print(x_c, y_c, z_c)
    return [x_c, y_c, z_c]


patch = 'inflow'
mesh_stats = readMeshStats()
face_labels = readFaceLabels(patch)
# print(face_labels)
face_points = readFacePoints(face_labels, mesh_stats)
point_labels = readPointLabels(mesh_stats['faces'])
point_cooordinates = processPointCoordinates(mesh_stats['points'])
point_centroids = processCentroids(face_labels, point_labels, point_cooordinates)
print(point_centroids)
# processCentroids(face_labels, point_labels, point_cooordinates)/
# print(face_points)