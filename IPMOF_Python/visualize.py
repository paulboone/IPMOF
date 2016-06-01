# Plots atom coordinates for packed unit cells
def plotPackedCell(packedCoor, azim, elev):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(packedCoor)):
        for coor in packedCoor[i]:
            r = i / len(packedCoor)
            ax.scatter(coor[0], coor[1], coor[2], '-o', c=[r, 0.1, 0.1], s=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()


# Plots atom coordinates for packed unit cells
def plotXYZ(xyzCoor, azim, elev):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for coor in xyzCoor:
        ax.scatter(coor[0], coor[1], coor[2], '-o', s=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()


def plotEnergyMap(eMap, azim, elev):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    xCoor = []
    yCoor = []
    zCoor = []
    energy = []
    for i in range(len(eMap)):
        xCoor.append(eMap[i][0])
        yCoor.append(eMap[i][1])
        zCoor.append(eMap[i][2])
        energy.append(eMap[i][3])
    scale = 1E20
    colors = []
    rgb = np.zeros([len(energy), 3])
    for i in range(len(energy)):
        if energy[i] > max(energy)*0.5/scale:
            energy[i] = 1
        elif energy[i] > max(energy)*0.25/scale:
            energy[i] = 0.9
        elif energy[i] > max(energy)*0.1/scale:
            energy[i] = 0.8
        elif energy[i] > max(energy)*0.075/scale:
            energy[i] = 0.7
        elif energy[i] > max(energy)*0.05/scale:
            energy[i] = 0.6
        elif energy[i] > max(energy)*0.025/scale:
            energy[i] = 0.5
        elif energy[i] > max(energy)*0.01/scale:
            energy[i] = 0.4
        else:
            energy[i] = 0.3
        colors.append(energy[i])
        rgb[i][0] = colors[i]

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(xCoor, yCoor, zCoor, 'o', c=rgb, edgecolor='')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()
