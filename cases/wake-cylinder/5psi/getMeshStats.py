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

    return mesh_stats 

print(readMeshStats())