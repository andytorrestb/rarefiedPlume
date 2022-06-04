from meshStats import readMeshStats
import numpy as np
import matplotlib.pyplot as plt

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

    # print(x_c, y_c, z_c)
    return [x_c, y_c, z_c]


def calculateNorm(face, point_labels, point_cooordinates):
    # Save needed cooridantes in a 2D list.
    face_vectors = []
    for i in range(len(point_labels[face])):
        coordinate = str(face_points[face][i])
        # print(coordinate, point_cooordinates[coordinate])
        position_vector = point_cooordinates[coordinate]
        face_vectors.append(position_vector)
        # print()
        # print()

    # Save relative poition vectors.
    Ax = float(face_vectors[0][0]) - float(face_vectors[1][0])
    Ay = float(face_vectors[0][1]) - float(face_vectors[1][1])
    Az = float(face_vectors[0][2]) - float(face_vectors[1][2])

    Bx = float(face_vectors[0][0]) - float(face_vectors[2][0])
    By = float(face_vectors[0][1]) - float(face_vectors[2][1])
    Bz = float(face_vectors[0][2]) - float(face_vectors[2][2])

    A = np.array([Ax, Ay, Az])
    B = np.array([Bx, By, Bz])
    return np.cross(A, B)

def processNorms(face_labels, point_labels, point_cooordinates):
    norm = {}

    for face in face_labels:
        face = str(face)
        norm[face] = calculateNorm(face, point_labels, point_cooordinates)
        print(norm[face])

    return norm

def graphNorm(face_norms, face_centroids):
    centroids = []
    norms = []

    # Convert data from dict to lists
    for face in face_centroids:
        centroids.append(face_centroids[str(face)])

    for face in face_norms:
        norms.append(face_norms[face])

    # Extract coordinate data, save to unique lists
    x_c = []
    y_c = []
    z_c = []
    for point in centroids:
        x_c.append(point[0])
        y_c.append(point[1])
        z_c.append(point[2])

    x_n = []
    y_n = []
    z_n = []
    for point in norms:
        x_n.append(point[0])
        y_n.append(point[1])
        z_n.append(point[2])

    # Plot data
    ax = plt.figure().add_subplot(projection = '3d')
    ax.quiver(x_c, y_c, z_c, x_n, y_n, z_n, length = 0.1, normalize = True)
    plt.show()
    # print(centroids)

patch = 'inflow'
mesh_stats = readMeshStats()
face_labels = readFaceLabels(patch)
face_points = readFacePoints(face_labels, mesh_stats)
point_labels = readPointLabels(mesh_stats['faces'])
point_cooordinates = processPointCoordinates(mesh_stats['points'])
face_centroids = processCentroids(face_labels, point_labels, point_cooordinates)
face_norms = processNorms(face_labels, point_labels, point_cooordinates)
graphNorm(face_norms, face_centroids)
