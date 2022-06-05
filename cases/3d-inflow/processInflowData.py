# from meshStats import readMeshStats
import numpy as np
import math
import matplotlib.pyplot as plt
import os

def readMeshStats():
    # Run check mesh, dump output inta log file. 
    os.system("checkMesh > log.checkMesh")

    mesh_stats = {}

    with open('log.checkMesh', 'r') as f:
        # print(f.read())
        line = f.readline()
        # Find mesh stats 
        while 'Mesh stats' not in line:
            line = f.readline()
        for i in range(9):
            line = f.readline()
            line = line.split(":")

            # Strip leading and trailing spaces. 
            for i in range(len(line)):
                line[i] = line[i].strip()

            # Replace internal space with underscores. 
            line[0] = line[0].replace(' ', '_')

            mesh_stats[line[0]] = int(line[1])

    os.system("rm log.checkMesh")

    return mesh_stats 

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

    with open(path, 'r') as f:
        # print(f.read())

        # Look for line containing nFaces
        line = f.readline()
        while str(nFaces) not in line:
            # print(str(nFaces))
            line = f.readline()

        # Skip a line
        line = f.readline()

        # Read the labeled IDs for each face in the patch.
        for i in range(nFaces):
            # face_labels.append(int(f.readline()))
            label = f.readline().split('\n')[0]
            face_labels[label] = int(label)
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

        # Use created dict to save labels for patch of interest to a seperate dict.
        for face in face_labels:
            # print(face)
            # print(face, face_dict[str(face)])
            # print(face_dict[str(face)].split("(")[1].split(")")[0].split(' '))
            face_points[str(face)] = face_dict[str(face)].split("(")[1].split(")")[0].split(' ')

    # Return faces for patch of interest.
    return face_points

def readPointLabels(tNFaces):
    path = 'constant/polyMesh/faces'

    # nFaces for the entire mesh.
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
    # Retrurn a dict that maps each label ID (str) to the correct coordinates.
    return point_labels

def processPointCoordinates(nPoints):
    # Read coordinates for all th points in the mesh.
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
    # Retrurn a dict that maps each label ID (str) to the correct coordinates.
    return point_cooordinates

def processCentroids(face_labels, point_labels, face_points, point_cooordinates):
    # Calculates centroids for each face in patch of interest.
    centroid = {}
    for face in face_labels:

        face = str(face)
        # print('face', face)
        # print('face_points[face]', face_points[face])
        # print('point_labels[face]', point_labels[face])
        # print(len(face_points[face]))

        # for i in range(len(face_points[face])):
        #     # print(point_cooordinates[str(face)])
        #     coordinate = str(face_points[face][i])
        #     print(coordinate, point_cooordinates[coordinate])

        # print()

        centroid[face] = calculateCentroid(face, face_points, point_cooordinates)
    # Return centroids as a dict that maps ID labels to the centroid vector.
    return centroid

def calculateCentroid(face, face_points, point_cooordinates):
    # Calculate and return the centroid vector for a given face
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


def calculateNorm(face, point_labels, face_points, point_cooordinates):
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
    return np.cross(B, A)

def processNorms(face_labels, point_labels, face_points, point_cooordinates):
    # Wrapper function for calcuteNorms
    norm = {}
    for face in face_labels:
        face = str(face)
        norm[face] = calculateNorm(face, point_labels, face_points, point_cooordinates)
        # print(norm[face])

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
    ax.set_xlim([0.0, 1.0])
    plt.show()
    # print(centroids)

def graphCentroid(face_centroids):
    centroids = []

    # Convert data from dict to lists
    for face in face_centroids:
        centroids.append(face_centroids[str(face)])

    # Extract coordinate data, save to unique lists
    x_c = []
    y_c = []
    z_c = []
    for point in centroids:
        x_c.append(point[0])
        y_c.append(point[1])
        z_c.append(point[2])

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x_c, y_c, z_c)
    ax.set_xlim([0.0, 1.0])
    plt.show()

def findSymmetryTheta(face_centroids):

    theta = {}

    for face in face_centroids:
        # print(type(face))
        print(face_centroids[face])
        y = float(face_centroids[face][1])
        z = float(face_centroids[face][2])
        theta[face] = math.atan(abs(y/z))

    return theta

def graphSymmetryTheta(face_centroids, face_thetas):
    centroids = []
    thetas = []

    # Convert data from dict to lists
    for face in face_centroids:
        centroids.append(face_centroids[str(face)])

    for face in face_thetas:
        thetas.append(face_thetas[str(face)])

    max_theta = max(thetas)

    for theta in thetas:
        print(theta)
        print()
        # theta = theta / max_theta
    # thetas = thetas / max(thetas)

    # Extract coordinate data, save to unique lists
    x_c = []
    y_c = []
    z_c = []
    for point in centroids:
        x_c.append(point[0])
        y_c.append(point[1])
        z_c.append(point[2])



    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x_c, y_c, z_c, c = thetas)
    ax.set_xlim([0.0, 1.0])
    plt.show()


def plumeSourceFlowModel():
    patch = 'inflow'
    mesh_stats = readMeshStats()
    face_labels = readFaceLabels(patch)
    face_points = readFacePoints(face_labels, mesh_stats)
    point_labels = readPointLabels(mesh_stats['faces'])
    point_cooordinates = processPointCoordinates(mesh_stats['points'])
    face_centroids = processCentroids(face_labels, point_labels, face_points, point_cooordinates)
    face_norms = processNorms(face_labels, point_labels, face_points, point_cooordinates)

    print('len(face_points)', len(face_points))
    print('len(face_labels)', len(face_labels))
    print('len(point_labels)', len(point_labels))
    print('len(point_cooordinates)', len(point_cooordinates))
    print('len(face_centroids)', len(face_centroids))
    print('len(face_norms)', len(face_norms))

    # graphNorm(face_norms, face_centroids)
    # graphCentroid(face_centroids

    face_thetas = findSymmetryTheta(face_centroids)
    graphSymmetryTheta(face_centroids, face_thetas)


plumeSourceFlowModel()