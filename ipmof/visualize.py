# IPMOF Visualization Functions
# Date: June 2016
# Author: Kutay B. Sezginel
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d


def plot_packed_cell(packed_uc_coors, edge_points, azim, elev):
    """
    Plots atom coordinates for packed unit cells.
     >>> plot_packed_cell(mof.packed_coors, mof.edge_points, 0, 0)
    """
    def orthogonal_proj(zfront, zback):
        a = (zfront + zback) / (zfront - zback)
        b = -2 * (zfront * zback) / (zfront - zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(packed_uc_coors)):
        for coor in packed_uc_coors[i]:
            r = i / len(packed_uc_coors)
            ax.scatter(coor[0], coor[1], coor[2], '-o', c=[r, 0.1, 0.1], s=10)

    for edge in edge_points:
        ax.scatter(edge[0], edge[1], edge[2], s=30)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()


def plot_unit_cell(atom_coors, edge_points, azim, elev):
    """
    Plots atom coordinates of a single unit cell with unit cell lines
     >>> plot_unit_cell(mof.atom_coors, mof.edge_points, 0, 0)

    view on a => z vs y (azim = 0 , elev = 180)
    view on c +> y vs x (axim = 90, elev = -90)
    """
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront + zback) / (zfront - zback)
        b = -2 * (zfront * zback) / (zfront - zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for coor in atom_coors:
        ax.scatter(coor[0], coor[1], coor[2], '-o', s=10)

    for edge in edge_points:
        ax.scatter(edge[0], edge[1], edge[2], c='r', s=30)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev
    # view on a => z vs y (azim = 0 , elev = 180)
    # view on b => x vs z => does not work (azim = 90 , elev = 0)
    # view on c +> y vs x (axim = 90, elev = -90)

    plt.show()


def plot_xyz(xyz_coor, azim, elev):
    """
    Plots atom coordinates of a single unit cell
     >>> plot_xyz(mof.atom_coors, 0, 0)
    """
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront + zback) / (zfront - zback)
        b = -2 * (zfront * zback) / (zfront - zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for coor in xyz_coor:
        ax.scatter(coor[0], coor[1], coor[2], '-o', s=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()


def plot_energy_map(emap, azim, elev):
    """
    Plots 3D color coded energy map.
     >>> plot_energy_map(emap, 0, 0)
    """
    def orthogonal_proj(zfront, zback):
        a = (zfront + zback) / (zfront - zback)
        b = -2 * (zfront * zback) / (zfront - zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])

    proj3d.persp_transformation = orthogonal_proj

    xCoor = []
    yCoor = []
    zCoor = []
    energy = []
    for i in range(len(emap)):
        xCoor.append(emap[i][0])
        yCoor.append(emap[i][1])
        zCoor.append(emap[i][2])
        energy.append(emap[i][3])
    scale = 1E20
    colors = []
    rgb = np.zeros([len(energy), 3])
    max_energy = max(energy)
    min_energy = min(energy)
    range_energy = abs(max_energy) + abs(min_energy)
    for i in range(len(energy)):
        if energy[i] > max(energy) * 0.5 / scale:
            energy[i] = 1
        elif energy[i] > max(energy) * 0.25 / scale:
            energy[i] = 0.9
        elif energy[i] > max(energy) * 0.1 / scale:
            energy[i] = 0.8
        elif energy[i] > max(energy) * 0.075 / scale:
            energy[i] = 0.7
        elif energy[i] > max(energy) * 0.05 / scale:
            energy[i] = 0.6
        elif energy[i] > max(energy) * 0.025 / scale:
            energy[i] = 0.5
        elif energy[i] > max(energy) * 0.01 / scale:
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
