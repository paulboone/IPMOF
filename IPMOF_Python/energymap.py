# IPMOF Energy Map Functions
# Date: June 2016
# Author: Kutay B. Sezginel
import os
import xlrd
import math
from math import floor
from math import ceil
from math import inf
import numpy as np
from forcefield import *


def energy_map(MOF1, atom_list, cut_off, grid_size):
    """
    Calculate energy map for a given MOF class.
    MOF -> base (map) | atomFFparameters -> sigma and epsilon values for given atoms
    cut_off -> cut-off value for LJ potential | grid_size -> grid size array for each dimension
    Packed coordinates for MOF1 must be defined before running the function.
    """
    sorted_x = sorted(MOF1.edgePoints, key=lambda x: x[0], reverse=True)
    sorted_y = sorted(MOF1.edgePoints, key=lambda y: y[1], reverse=True)
    sorted_z = sorted(MOF1.edgePoints, key=lambda z: z[2], reverse=True)
    emap_max = [ceil(sorted_x[0][0]), ceil(sorted_y[0][1]), ceil(sorted_z[0][2])]
    emap_min = [floor(sorted_x[-1][0]), floor(sorted_y[-1][1]), floor(sorted_z[-1][2])]

    x_grid = np.linspace(emap_min[0], emap_max[0], (emap_max[0] - emap_min[0]) / grid_size + 1)
    y_grid = np.linspace(emap_min[1], emap_max[1], (emap_max[1] - emap_min[1]) / grid_size + 1)
    z_grid = np.linspace(emap_min[2], emap_max[2], (emap_max[2] - emap_min[2]) / grid_size + 1)

    num_atoms = len(atom_list['sigma'])

    # Initialize energy map according to grid size and coordinates plus number of unique atoms
    energy_map = np.zeros([len(x_grid) * len(y_grid) * len(z_grid), num_atoms + 3])

    sig, eps = LBmix(MOF1.sigma, atom_list['sigma'], MOF1.epsilon, atom_list['epsilon'])

    map_index = 0
    v = np.zeros([num_atoms])

    for x in x_grid:
        for y in y_grid:
            for z in z_grid:
                energy_map[map_index][0:3] = [x, y, z]
                v_total = np.zeros([num_atoms])
                for unit_cell in MOF1.packedCoor:
                    MOF1index = 0
                    for atomCoor in unit_cell:
                        atom_index_1 = MOF1.uniqueAtomNames.index(MOF1.atomName[MOF1index])
                        MOF1index += 1
                        dist = coor_dist(atomCoor, [x, y, z])
                        if dist > cut_off:
                            continue
                        if dist == 0:
                            energy_map[map_index][3:(num_atoms + 3)] = np.ones([1, num_atoms]) * inf
                        else:
                            for atom_index_2 in range(num_atoms):
                                sig_mix = sig[atom_index_1][atom_index_2]
                                eps_mix = eps[atom_index_1][atom_index_2]
                                v[atom_index_2] = lennard_jones(dist, sig_mix, eps_mix)
                                v_total[atom_index_2] = v_total[atom_index_2] + v[atom_index_2]
                energy_map[map_index][3:(num_atoms + 3)] = v_total
                map_index += 1
    return energy_map


def coor_dist(coor1, coor2):
    """
    Calculates distance between two given coordinates: [x1, y1, z1] and [x2, y2, z2]
    """
    return math.sqrt((coor1[0] - coor2[0])**2 + (coor1[1] - coor2[1])**2 + (coor1[2] - coor2[2])**2)
