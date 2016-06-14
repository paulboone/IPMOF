# IPMOF Energy Map Functions
# Date: June 2016
# Author: Kutay B. Sezginel
import os
from math import floor, ceil, inf, sqrt

import xlrd
import numpy as np

from ipmof.forcefield import lorentz_berthelot_mix, lennard_jones
from ipmof.crystal import MOF


def energy_map(sim_par, mof, atom_list):
    """
    Calculate energy map for a given MOF class with following properties:
        - edge_points   - uniq_atom_names   - atom_names    - packed_coors
        - sigma         - epsilon
    MOF -> base (map) | atomFFparameters -> sigma and epsilon values for given atoms
    cut_off -> cut-off value for LJ potential | grid_size -> grid size array for each dimension
    Packed coordinates for MOF must be defined before running the function.
    Resulting energy map is structured as follows:
        emap[0] = [x, y, z, atom1_energy, atom2_energy, atom3_energy, ...]
    """
    cut_off = sim_par['cut_off']
    grid_size = sim_par['grid_size']
    # Determine max and min coordinates for the unit cell to construct bounding box grid
    sorted_x = sorted(mof.edge_points, key=lambda x: x[0], reverse=True)
    sorted_y = sorted(mof.edge_points, key=lambda y: y[1], reverse=True)
    sorted_z = sorted(mof.edge_points, key=lambda z: z[2], reverse=True)
    emap_max = [ceil(sorted_x[0][0]), ceil(sorted_y[0][1]), ceil(sorted_z[0][2])]
    emap_min = [floor(sorted_x[-1][0]), floor(sorted_y[-1][1]), floor(sorted_z[-1][2])]

    x_grid = np.linspace(emap_min[0], emap_max[0], (emap_max[0] - emap_min[0]) / grid_size + 1)
    y_grid = np.linspace(emap_min[1], emap_max[1], (emap_max[1] - emap_min[1]) / grid_size + 1)
    z_grid = np.linspace(emap_min[2], emap_max[2], (emap_max[2] - emap_min[2]) / grid_size + 1)

    num_atoms = len(atom_list['sigma'])

    # Initialize energy map according to grid size and coordinates plus number of unique atoms
    energy_map = np.zeros([len(x_grid) * len(y_grid) * len(z_grid), num_atoms + 3])

    sig, eps = lorentz_berthelot_mix(mof.sigma, atom_list['sigma'], mof.epsilon, atom_list['epsilon'])

    map_index = 0
    v = np.zeros([num_atoms])

    for x in x_grid:
        for y in y_grid:
            for z in z_grid:
                energy_map[map_index][0:3] = [x, y, z]
                v_total = np.zeros([num_atoms])
                for unit_cell in mof.packed_coors:
                    mof_index = 0
                    for atom_coor in unit_cell:
                        atom_index_1 = mof.uniq_atom_names.index(mof.atom_names[mof_index])
                        mof_index += 1
                        dist = coor_dist(atom_coor, [x, y, z])
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


def energy_map_index(coor, emap_max, emap_min):
    """
    Finds index of a coordinate in the energy map. Only works for grid_size = 1.
    emap_max = [emap[-1][0], emap[-1][1], emap[-1][2]]
    emap_min = [emap[0][0], emap[0][1], emap[0][2]]
    """
    side_length = []
    round_coor = []

    for c, emax, emin in zip(coor, emap_max, emap_min):
        side_length.append(emax - emin + 1)
        round_coor.append(round(c))
    emap_index = round_coor[0] * side_length[1] * side_length[2]
    emap_index += round_coor[1] * side_length[2] + round_coor[2]

    return int(emap_index)


def energy_map_atom_index(atom_name, emap_atom_list):
    """
    Returns index of a given atom in the energy map.

    given an energy map in the form:
        emap[i] = [x, y, z, C_atom_energy, H_atom_energy, O_atom_energy, Zn_atom_energy]

    Example usage:
     >>> atom_index = energy_map_atom_index('O', atom_list)
     ... 5

    so emap[i][5] would give the energy value for O atom
    """
    for atom_index, atom in enumerate(emap_atom_list['atom']):
        if atom == atom_name:
            emap_atom_index = atom_index

    return int(emap_atom_index + 3)


def coor_dist(coor1, coor2):
    """
    Calculates distance between two given coordinates: [x1, y1, z1] and [x2, y2, z2]
    """
    return sqrt((coor1[0] - coor2[0])**2 + (coor1[1] - coor2[1])**2 + (coor1[2] - coor2[2])**2)


def get_mof_file_list(file_dir, file_format):
    """
    Generates a list of MOF file names in a given directory and MOF file format
    """
    file_list = os.listdir(file_dir)
    mof_list = []
    for file_name in file_list:
        if file_format in file_name:
            mof_list.append(file_name)

    return mof_list


def get_mof_list(mof_path_list, force_field):
    """
    Generates a list of MOF objects using a given list of MOF file directories
    """
    mof_list = []
    for mof in mof_path_list:
        mof_obj = MOF(mof)
        mof_obj.force_field(force_field)
        mof_list.append(mof_obj)

    return mof_list


def get_uniq_atom_list(mof_list):
    """
    Gets atom name, epsilon, and sigma values for non-repeating (unique) atoms in a list of
    MOF classes.
    """
    all_atom_list = {'atom': [], 'sigma': [], 'epsilon': []}
    for mof in mof_list:
        for atom, sig, eps in zip(mof.uniq_atom_names, mof.sigma, mof.epsilon):
            all_atom_list['atom'].append(atom)
            all_atom_list['sigma'].append(sig)
            all_atom_list['epsilon'].append(eps)

    uniq_atom_list = {'atom': [], 'sigma': [], 'epsilon': []}
    # Creates a list of unique entries in all_atom_list atom names
    uniq_atom_list['atom'] = list(set(all_atom_list['atom']))
    # Initializes epsilon and sigma arrays according to the length of unique atom names
    uniq_atom_list['epsilon'] = [0] * len(uniq_atom_list['atom'])
    uniq_atom_list['sigma'] = [0] * len(uniq_atom_list['atom'])

    # Finds epsilon and sigma values for unique atoms by comparing the atom names in both lists
    # Epsilon and sigma values are overwritten multiple times for repeating atoms
    for atom, sig, eps in zip(all_atom_list['atom'], all_atom_list['sigma'], all_atom_list['epsilon']):
        if atom in uniq_atom_list['atom']:
            uniq_index = uniq_atom_list['atom'].index(atom)
            uniq_atom_list['epsilon'][uniq_index] = eps
            uniq_atom_list['sigma'][uniq_index] = sig

    return uniq_atom_list
