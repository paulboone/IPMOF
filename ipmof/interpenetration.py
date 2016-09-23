# IPMOF interpenetration library
# Date: June 2016
# Author: Kutay B. Sezginel
import os
import math
from random import random
import shutil
from glob import glob

from ipmof.crystal import Packing, MOF
from ipmof.geometry import Coor, Quaternion
from ipmof.energymap import energy_map_index, energy_map_atom_index, import_energy_map, get_mof_list
from ipmof.parameters import export_interpenetration_results


def initial_coordinates(mof, energy_map, atom_list, energy_limit):
    """
    Determine initial coordinates to start interpenetration simulations from.
    Points are determined acoording to their energy (accepted if energy < energy_limit)
    and their position (accepted if applying pbc does not change its coordinates)
    """
    reference_atom = 'C'
    if reference_atom in atom_list['atom']:
        ref_atom_index = int(atom_list['atom'].index(reference_atom) + 3)
    else:
        ref_atom_index = 3
    initial_coors = []
    energy_count = 0
    pbc_count = 0
    for emap_line in energy_map:
        emap_coor = Coor([emap_line[0], emap_line[1], emap_line[2]])
        pbc_coor = emap_coor.pbc(mof.uc_size, mof.uc_angle, mof.frac_ucv)
        pbc_x = round(pbc_coor.x, 1)
        pbc_y = round(pbc_coor.y, 1)
        pbc_z = round(pbc_coor.z, 1)
        # print(emap_coor.x, pbc_x)
        if pbc_x == emap_coor.x and pbc_y == emap_coor.y and pbc_z == emap_coor.z:
            if emap_line[ref_atom_index] < energy_limit:
                initial_coors.append(Coor([emap_line[0], emap_line[1], emap_line[2]]))
            else:
                energy_count += 1
        else:
            pbc_count += 1

    # print('Ommited PBC: ', pbc_count, ' Energy: ', energy_count)
    return initial_coors


def tripolate(point, atom_index, emap, x_length, y_length):
    """
    3D Linear Interpolation for given energy map and point in space (point must be in emap).
    Only works for energy map constructed with a grid size of 1.
    """
    point0 = []
    dif = []
    for p in point:
        point0.append(math.floor(p))
        dif.append(p - point0[-1])

    i000 = int(point0[0] * x_length + point0[1] * y_length + point0[2])
    i001 = i000 + 1
    i010 = i000 + y_length
    i011 = i010 + 1
    i100 = i000 + x_length
    i101 = i100 + 1
    i110 = i010 + x_length
    i111 = i110 + 1

    d1 = 1 - dif[0]
    c00 = emap[i000][atom_index] * d1 + emap[i100][atom_index] * dif[0]
    c01 = emap[i001][atom_index] * d1 + emap[i101][atom_index] * dif[0]
    c10 = emap[i010][atom_index] * d1 + emap[i110][atom_index] * dif[0]
    c11 = emap[i011][atom_index] * d1 + emap[i111][atom_index] * dif[0]

    c0 = c00 * (1 - dif[1]) + c10 * dif[1]
    c1 = c01 * (1 - dif[1]) + c11 * dif[1]

    c = c0 * (1 - dif[2]) + c1 * dif[2]

    return c


def check_interpenetration(sim_par, base_mof, mobile_mof, emap, atom_list):
    """
    Run interpenetration algorithm with given simulation parameters and energy map.
    Returns simulation summary and structural information on the discovered structures.
    """
    # Initialize simulation parameters
    structure_energy_limit = sim_par['structure_energy_limit']
    atom_energy_limit = sim_par['atom_energy_limit']
    rotation_limit = sim_par['rotation_limit']
    rotation_freedom = sim_par['rotation_freedom']
    summary_percent = sim_par['summary_percent']

    emap_max = [emap[-1][0], emap[-1][1], emap[-1][2]]
    emap_min = [emap[0][0], emap[0][1], emap[0][2]]
    side_length = [emap_max[0] - emap_min[0] + 1, emap_max[1] - emap_min[1] + 1, emap_max[2] - emap_min[2] + 1]
    x_length, y_length = int(side_length[1] * side_length[2]), int(side_length[2])

    # Initialize simulation variables
    Quat = Quaternion([0, 1, 1, 1])

    initial_coors = initial_coordinates(base_mof, emap, atom_list, atom_energy_limit)
    trial_limit = len(initial_coors) * rotation_limit
    rot_freedom = 360 / rotation_freedom
    # omitted_coordinates = len(emap) - len(initial_coors)

    div = round(trial_limit / (100 / summary_percent))
    summary = {'percent': [], 'structure_count': [], 'trial_count': []}
    new_structures = []

    abort_ip = False
    structure_count = 0
    structure_total_energy = 0
    initial_coor_index = 0

    # Main interpenetration algorithm
    for t in range(trial_limit):  # Can iterate over something else???
        abort_ip = False
        # Interpenetration trial loop
        # Try interpenetration for a specific orientation by going through each atom in mobile mof
        for idx in range(len(mobile_mof)):    # Can iterate over something else???

            if not abort_ip:
                # If the interpenetration is just starting select rotation angles
                if idx == 0:
                    if t % rotation_limit == 0:
                        first_point = initial_coors[initial_coor_index]
                        initial_coor_index += 1
                        # initial_coor_trial_count += 1

                    # Determine random angles for rotation in 3D space
                    x_angle = 2 * math.pi * math.floor(random() * rot_freedom) / rot_freedom
                    y_angle = 2 * math.pi * math.floor(random() * rot_freedom) / rot_freedom
                    z_angle = 2 * math.pi * math.floor(random() * rot_freedom) / rot_freedom

                    # Rotate first atom of the mobile MOF
                    atom_name = mobile_mof.atom_names[idx]
                    new_coor = Coor(mobile_mof.atom_coors[idx])
                    Q = Quaternion([1, new_coor.x, new_coor.y, new_coor.z])  # Might be a better way to do this
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [1, 0, 0], x_angle)
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 1, 0], y_angle)
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 0, 1], z_angle)
                    # new_coor = Q.coor()
                    new_coor = Coor(Q.xyz())

                    translation_vector = first_point - new_coor  # Check if operation is correct

                    # Initialize new structure dictionay
                    structure = {'atom_names': [], 'atom_coors': [],
                                 'pbc_coors': [], 'energy': [], 'rotation': []}
                    structure['first_point'] = first_point.xyz()
                    structure['translation_vector'] = translation_vector.xyz()
                    structure['atom_coors'].append(first_point.xyz())  # Why first point not new_coor?
                    structure['pbc_coors'].append(first_point.xyz())
                    structure['atom_names'].append(atom_name)
                    structure['rotation'] = [x_angle, y_angle, z_angle]

                # If interpenetration is still going on
                if idx < len(base_mof) and idx > 0:
                    atom_name = atom_name = mobile_mof.atom_names[idx]
                    new_coor = Coor(mobile_mof.atom_coors[idx])
                    Q = Quaternion([1, new_coor.x, new_coor.y, new_coor.z])  # Might be a better way to do this
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [1, 0, 0], x_angle)
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 1, 0], y_angle)
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 0, 1], z_angle)
                    # new_coor = Q.coor()
                    new_coor = Coor(Q.xyz())

                    new_coor += translation_vector
                    pbc_coor = new_coor.pbc(base_mof.uc_size, base_mof.uc_angle, base_mof.frac_ucv)

                    emap_index = energy_map_index(pbc_coor.xyz(), emap_max, emap_min)
                    emap_atom_index = energy_map_atom_index(atom_name, atom_list)

                    point_energy = tripolate(pbc_coor.xyz(), emap_atom_index, emap, x_length, y_length)
                    structure_total_energy += point_energy

                    if structure_total_energy > structure_energy_limit:
                        structure_total_energy = 0
                        abort_ip = True
                        break  # Fix this part (break interpenetration trial loop)
                    elif point_energy > atom_energy_limit:
                        structure_total_energy = 0
                        abort_ip = True
                        break  # Fix this part (break interpenetration trial loop)
                    else:
                        structure['atom_coors'].append(new_coor.xyz())
                        structure['pbc_coors'].append(pbc_coor.xyz())
                        structure['atom_names'].append(atom_name)

                # If interpenetration trial ended with no collision
                if idx == len(mobile_mof) - 1:
                    # Record structure information!!!!
                    structure['energy'] = structure_total_energy
                    new_structures.append(structure)
                    structure_count += 1
                    structure_total_energy = 0

        # Record simulation progress according to division (div)
        if t % div == 0:
            percent_complete = round(t / trial_limit * 100)

            # Record summary information
            summary['percent'].append(percent_complete)
            summary['structure_count'].append(structure_count)
            summary['trial_count'].append(t)

    return summary, new_structures


def check_extension(sim_par, base_mof, mobile_mof, emap, emap_atom_list, new_structure):
    """
    *** Not Complete ***
    Checks collision between interpenetrating layer and base layer for a determined distance.
    Distance is calculated from given ext_cut_off value which determines the packing amount of the
    interpenetrating layer.
    Each coordinate in the interpenetrating layer is checked for high energy values by applying
    perodic boundary conditions to the coordinate according to energy map of the base layer.
    """
    emap_max = [emap[-1][0], emap[-1][1], emap[-1][2]]
    emap_min = [emap[0][0], emap[0][1], emap[0][2]]
    side_length = [emap_max[0] - emap_min[0] + 1, emap_max[1] - emap_min[1] + 1, emap_max[2] - emap_min[2] + 1]
    x_length, y_length = int(side_length[1] * side_length[2]), int(side_length[2])

    energy_limit = sim_par['atom_energy_limit']
    ext_cut_off = sim_par['ext_cut_off']

    rotation_info = new_structure['rotation']
    first_point = new_structure['first_point']
    translation_vector = Coor(new_structure['translation_vector'])

    Quat = Quaternion([0, 1, 1, 1])

    packing_factor = Packing.factor(mobile_mof.uc_size, ext_cut_off)
    uc_vectors = Packing.uc_vectors(mobile_mof.uc_size, mobile_mof.uc_angle)
    trans_vec = Packing.translation_vectors(packing_factor, uc_vectors)
    packed_coors = Packing.uc_coors(trans_vec, packing_factor, uc_vectors, mobile_mof.atom_coors)

    x_angle = rotation_info[0]
    y_angle = rotation_info[1]
    z_angle = rotation_info[2]

    collision = False
    collision_info = {'exist': collision, 'coor': None, 'pbc_coor': None}
    for unit_cell in packed_coors:

        if not collision:

            for coor_index, coor in enumerate(unit_cell):

                if not collision:

                    new_coor = Coor(coor)
                    Q = Quaternion([1, new_coor.x, new_coor.y, new_coor.z])  # Might be a better way to do this
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [1, 0, 0], x_angle)
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 1, 0], y_angle)
                    Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 0, 1], z_angle)
                    new_coor = Coor(Q.xyz())

                    new_coor += translation_vector
                    pbc_coor = new_coor.pbc(base_mof.uc_size, base_mof.uc_angle, base_mof.frac_ucv)

                    atom_name = mobile_mof.atom_names[coor_index]
                    atom_index = energy_map_atom_index(atom_name, emap_atom_list)

                    point_energy = tripolate(pbc_coor.xyz(), atom_index, emap, x_length, y_length)

                    if point_energy < energy_limit:
                        continue
                    else:
                        collision = True
                        collision_info = {'exist': collision, 'coor': new_coor.xyz(), 'pbc_coor': pbc_coor.xyz()}
                        break
                else:
                    break
        else:
            break

    return collision_info


def save_extension(sim_par, base_mof, mobile_mof, emap, emap_atom_list, new_structure):
    """
    Using the rotation_info and translation_vector from interpenetration to extended_coors
    mobile structure for a given distance (cut_off).
    Returns atom names, coordinates and packing factor.
    """
    emap_max = [emap[-1][0], emap[-1][1], emap[-1][2]]
    emap_min = [emap[0][0], emap[0][1], emap[0][2]]

    export_cut_off = sim_par['cut_off']

    rotation_info = new_structure['rotation']
    first_point = new_structure['first_point']
    translation_vector = Coor(new_structure['translation_vector'])

    Quat = Quaternion([0, 1, 1, 1])

    packing_factor = Packing.factor(mobile_mof.uc_size, export_cut_off)
    uc_vectors = Packing.uc_vectors(mobile_mof.uc_size, mobile_mof.uc_angle)
    trans_vec = Packing.translation_vectors(packing_factor, uc_vectors)
    packed_coors = Packing.uc_coors(trans_vec, packing_factor, uc_vectors, mobile_mof.atom_coors)

    x_angle = rotation_info[0]
    y_angle = rotation_info[1]
    z_angle = rotation_info[2]

    extended_coors = []
    extended_names = []

    for unit_cell in packed_coors:

        for coor_index, coor in enumerate(unit_cell):

            new_coor = Coor(coor)
            Q = Quaternion([1, new_coor.x, new_coor.y, new_coor.z])  # Might be a better way to do this
            Q = Quat.rotation(Q.xyz(), [0, 0, 0], [1, 0, 0], x_angle)
            Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 1, 0], y_angle)
            Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 0, 1], z_angle)
            new_coor = Coor(Q.xyz())

            new_coor += translation_vector

            atom_name = mobile_mof.atom_names[coor_index]

            extended_names.append(atom_name)
            extended_coors.append(new_coor.xyz())

    extended_structure = {'atom_names': extended_names, 'atom_coors': extended_coors}
    extended_structure['name'] = mobile_mof.name
    extended_structure['packing_factor'] = packing_factor

    return extended_structure


def run_interpenetration(interpenetration_path, sim_par, sim_dir):
    """
    Interpenetration algorithm for job server.
    1) Checks Interpenetration
    2) Gets minimum energy structures
        - Performs collision check by extending interpenetrating structure
        - Saves requested structure files
    """
    emap_path, base_mof_path, mobile_mof_path = interpenetration_path
    base_mof = MOF(base_mof_path)
    mobile_mof = MOF(mobile_mof_path)
    atom_list, emap = import_energy_map(emap_path)

    summary, new_structures = check_interpenetration(sim_par, base_mof, mobile_mof, emap, atom_list)

    # Create export directory ------------------=-------------------------------------------
    export_dir = os.path.join(sim_dir['export_dir'], base_mof.name + '_' + mobile_mof.name)
    if os.path.exists(export_dir):
        shutil.rmtree(export_dir)
    os.mkdir(export_dir)

    structure_info = [{'S1': base_mof.name, 'S2': mobile_mof.name, 'Structures': len(new_structures)}]
    # Export Min Energy Structures ---------------------------------------------------------
    if len(new_structures) > 0:

        export_count = min(len(new_structures), sim_par['export_structures'])
        for export_index in range(export_count):
            # Get minimum energy structure by sorting total structure energies
            min_energy_structure = sorted(new_structures, key=lambda k: k['energy'])[export_index]

            # Check for collision in the extended unitcell of new structure and energy map
            collision = check_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)

            # Record structure information -------------------------------------------------
            # Float is used to convert values to numbers for proper storage with yaml format
            structure_info.append({'energy': float(round(min_energy_structure['energy'], 3)),
                                   'collision': collision,
                                   'rotation': [round(math.degrees(a)) for a in min_energy_structure['rotation']],
                                   'initial_coordinate': [float(round(p, 1)) for p in min_energy_structure['first_point']]})

            # Export new structure(s) -------------------------------------------------------
            export_structures(sim_par, base_mof, mobile_mof, min_energy_structure, emap, atom_list, export_index, export_dir)

    export_interpenetration_results(sim_par, structure_info, summary, export_dir)


def get_interpenetration_list(sim_par, sim_dir):
    """
    Returns a dictionary containing paths for energy map and mof files for each MOF combination.
    Energy map is read from ~/energy_map and MOFs are read from ~/mof.
    All possible combinations are then arranged in a dictionary format:
        - If self_interpenetration parameter if False homo-interpenetration is excluded
        - Reverse combinations are excluded (ex: if MOF1_MOF2 in list then MOF2_MOF1 is not selected)

    Format: interpenetration_list = {'emap_path': [], 'emap_mof_path': [], 'ip_mof_path': []}
    """
    mof_path_list = get_mof_list(sim_par, sim_dir)
    emap_path_list = os.listdir(sim_dir['energy_map_dir'])

    interpenetration_list = []
    ip_mof_list = []
    emap_mof_list = []

    for emap_path in emap_path_list:

        emap_mof_name = os.path.basename(emap_path).split('_emap')[0]
        # What if multiple files are returned?? (IndexError if file cannot be found)
        emap_mof_path = glob(os.path.join(sim_dir['mof_dir'], emap_mof_name) + '*')[0]
        emap_path = os.path.join(sim_dir['energy_map_dir'], emap_path)

        for ip_mof_path in mof_path_list:

            if emap_mof_path == ip_mof_path and not sim_par['self_interpenetration']:
                continue

            elif emap_mof_path in ip_mof_list and ip_mof_path in emap_mof_list and emap_mof_path != ip_mof_path:
                continue

            else:
                interpenetration_list.append((emap_path, emap_mof_path, ip_mof_path))
                ip_mof_list.append(ip_mof_path)
                emap_mof_list.append(emap_mof_path)

    return interpenetration_list


def export_structures(sim_par, base_mof, mobile_mof, min_energy_structure, emap, atom_list, export_index, export_dir):
    """
    Export requested interpenetration structures.
    Types of structures can be selected from simulation parameters.
    """
    new_structure = {'atom_names': min_energy_structure['atom_names'], 'name': mobile_mof.name}
    if sim_par['export_pbc']:
        new_structure['atom_coors'] = min_energy_structure['pbc_coors']
    else:
        new_structure['atom_coors'] = min_energy_structure['atom_coors']

    if sim_par['export_single']:
        new_mobile_mof = MOF(new_structure, file_format='dict')
        joined_mof = base_mof.join(new_mobile_mof, colorify=False)
        joined_mof.name += '_' + str(export_index + 1)
        joined_mof.export(export_dir, file_format=sim_par['export_format'])

    if sim_par['export_single_color']:
        # Join base and mobile structure layers
        new_mobile_mof = MOF(new_structure, file_format='dict')
        joined_mof_color = base_mof.join(new_mobile_mof, colorify=True)
        joined_mof_color.name += '_' + str(export_index + 1) + 'C'
        joined_mof_color.export(export_dir, file_format=sim_par['export_format'])

    if sim_par['export_packed']:
        # Pack new structure by using rotation and first point information
        extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])
        packed_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
        joined_structure = {'atom_coors': packed_structure['atom_coors'] + extended_structure['atom_coors'],
                            'atom_names': packed_structure['atom_names'] + extended_structure['atom_names'],
                            'name': packed_structure['name'] + '_' + extended_structure['name']}
        joined_packed_mof = MOF(joined_structure, file_format='dict')
        joined_packed_mof.name += '_' + str(export_index + 1) + 'P'
        joined_packed_mof.export(export_dir, file_format='xyz')

    if sim_par['export_packed_color']:
        # Pack new structure by using rotation and first point information
        extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])
        packed_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
        joined_structure = {'atom_coors': packed_structure['atom_coors'] + extended_structure['atom_coors'],
                            'atom_names': ['C'] * len(packed_structure['atom_names']) +
                                          ['O'] * len(extended_structure['atom_names']),
                            'name': packed_structure['name'] + '_' + extended_structure['name']}
        joined_packed_mof = MOF(joined_structure, file_format='dict')
        joined_packed_mof.name += '_' + str(export_index + 1) + 'PC'
        joined_packed_mof.export(export_dir, file_format='xyz')
