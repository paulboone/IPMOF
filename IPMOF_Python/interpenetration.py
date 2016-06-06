# Quaternion class
# - Multiplication  - Division  - Inverse   - Rotation  - Output List <Q.xyz()>
#  >>> Q = Quaternion([0, 1, 1, 1])
# Date: June 2016
# Author: Kutay B. Sezginel
import math
from crystal import Packing, MOF
from quaternion import Quaternion
from energymap import energy_map_index, energy_map_atom_index


class Coor(object):
    """
    Coor class for holding 3D space coordinates.
    """
    def __init__(self, input):
        if isinstance(input, list):
            self.x = input[0]
            self.y = input[1]
            self.z = input[2]
        else:
            raise TypeError('Input type not supported. Use list [x, y, z]')

    def __repr__(self):
        return "<Coordinate object x:%s y:%s z:%s>" % (self.x, self.y, self.z)

    def __str__(self):
        return "%s %s %s" % (round(self.x, 4), round(self.y, 4), round(self.z, 4))

    def __add__(self, coor2):
        return Coor([self.x + coor2.x, self.y + coor2.y, self.z + coor2.z])

    def __sub__(self, coor2):
        return Coor([self.x - coor2.x, self.y - coor2.y, self.z - coor2.z])

    def dist(self, coor2):
        """
        Calculates distance between this coordinate and another given coordinate.

        Example usage:
         >>> coor1 = Coor([1, 2, 3])
         >>> coor2 = Coor([-3, 4, -2])
         >>> coor1.dist(coor2) -> 6.708203932499369
        """
        dx = self.x - coor2.x
        dy = self.y - coor2.y
        dz = self.z - coor2.z

        return math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    def xyz(self):
        """
        Return a list containing x, y, z coordinates.

        Example usage:
         >>> coor1 = Coor([1, 2, 3])
         >>> coor1.xyz() -> [1, 2, 3]
        """
        return [self.x, self.y, self.z]

    def frac(self, uc_size, uc_angle, frac_ucv):
        """
        Converts cartesian coordinates to fractional coordinates.
        *The fractional unit cell volume is calculated each time. Instead can be given as input.

        Example usage:
         >>> coor1 = Coor([-11, 22, 33])
         >>> coor1.frac([26, 26, 26], [90, 90, 90]) -> <Coordinate object x:-0.423 y:0.846 z:1.269>
        """
        alp = math.radians(uc_angle[0])
        bet = math.radians(uc_angle[1])
        gam = math.radians(uc_angle[2])

        v = frac_ucv
        # v = 1 - math.cos(alp) ** 2 - math.cos(bet) ** 2
        # v += - math.cos(gam) ** 2 + 2 * math.cos(alp) * math.cos(bet) * math.cos(gam)
        # v = math.sqrt(v)

        a = uc_size[0]
        b = uc_size[1]
        c = uc_size[2]

        x = self.x
        y = self.y
        z = self.z

        x_frac = 1 / a * x
        x_frac += - math.cos(gam) / (a * math.sin(gam)) * y
        x_frac += (math.cos(alp) * math.cos(gam) - math.cos(bet)) / (a * v * math.sin(gam)) * z

        y_frac = 1 / (b * math.sin(gam)) * y
        y_frac += (math.cos(bet) * math.cos(gam) - math.cos(alp)) / (b * v * math.sin(gam)) * z

        z_frac = math.sin(gam) / (c * v) * z

        return Coor([x_frac, y_frac, z_frac])

    def car(self, uc_size, uc_angle, frac_ucv):
        """
        Converts fractional coordinates to cartesian coordinates.
        Takes unit cell size, angles and fractional unit cell volume as input.

        Example usage:
         >>> coor1 = Coor([0.23, 2.21, -1.33])
         >>> coor1.car([26, 26, 26], [90, 90, 90]) -> <Coordinate object x:5.98 y:57.46 z:-34.58>
        """
        alp = math.radians(uc_angle[0])
        bet = math.radians(uc_angle[1])
        gam = math.radians(uc_angle[2])

        v = frac_ucv
        # v = 1 - math.cos(alp) ** 2 - math.cos(bet) ** 2
        # v += - math.cos(gam) ** 2 + 2 * math.cos(alp) * math.cos(bet) * math.cos(gam)
        # v = math.sqrt(v)

        a = uc_size[0]
        b = uc_size[1]
        c = uc_size[2]

        x_frac = self.x
        y_frac = self.y
        z_frac = self.z

        x = a * x_frac
        x += b * math.cos(gam) * y_frac
        x += c * math.cos(bet) * z_frac

        y = b * math.sin(gam) * y_frac
        y += c * (math.cos(alp) - math.cos(bet) * math.cos(gam)) / math.sin(gam) * z_frac

        z = c * v / math.sin(gam) * z_frac

        return Coor([x, y, z])

    def pbc(self, uc_size, uc_angle, frac_ucv):
        """
        Apply perodic boundary conditions to given cartesian coordinates and unit cell parameters.

        Example usage:
         >>> coor1 = Coor([0.23, 2.21, -1.33])
         >>> coor1.pbc([26, 26, 26], [90, 90, 90]) -> <Coordinate object x:0.23 y:0.21 z:0.67>
        """
        frac_coor = self.frac(uc_size, uc_angle, frac_ucv)
        frac_pbc_coor = frac_coor.frac_pbc()
        car_pbc_coor = frac_pbc_coor.car(uc_size, uc_angle, frac_ucv)

        return car_pbc_coor

    def frac_pbc(self):
        """
        Apply perodic boundary conditions to given fractional coordinates.

        Example usage:
         >>> coor1 = Coor([0.23, 2.21, -1.33])
         >>> coor1.frac_pbc() -> <Coordinate object x:0.23 y:0.21 z:0.67>
        """
        pbc_x = self.x - math.floor(self.x)
        pbc_y = self.y - math.floor(self.y)
        pbc_z = self.z - math.floor(self.z)

        return Coor([pbc_x, pbc_y, pbc_z])


def initial_coordinates(MOF, energy_map, atom_list, energy_limit):
    """
    Determine initial coordinates to start interpenetration simulations from.
    Points are determined acoording to their energy (accepted if energy < energy_limit)
    and their position (accepted if applying pbc does not change its coordinates)
    """
    reference_atom = 'C'
    ref_atom_index = atom_list['atom'].index(reference_atom) + 3
    initial_coors = []
    energy_count = 0
    pbc_count = 0
    for emap_line in energy_map:
        emap_coor = Coor([emap_line[0], emap_line[1], emap_line[2]])
        pbc_coor = emap_coor.pbc(MOF.uc_size, MOF.uc_angle, MOF.frac_ucv)
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


def trilinear_interpolate(point, atom_index, emap, emap_max, emap_min):
    """
    *** Working but needs to be checked for accuracy ***
    3D Linear Interpolation for given energy map and point in space (point must be in emap).
    Only works for energy map constructed with a grid size of 1.
    MODIFIES INPUT COORDINATE IF A ROUNDED COORDINATE VALUE IS GIVEN!!!!!!!!!!!!!
    """
    point1 = []
    point0 = []
    dif = []
    for p in point:
        if round(p) == p:
            p += 1E-10
        point0.append(math.floor(p))
        point1.append(math.ceil(p))
        dif.append((p - point0[-1]) / (point1[-1] - point0[-1]))

    i000 = energy_map_index(point0, emap_max, emap_min)                               # (0, 0, 0)
    i100 = energy_map_index([point1[0], point0[1], point0[2]], emap_max, emap_min)    # (1, 0, 0)
    i001 = energy_map_index([point0[0], point0[1], point1[2]], emap_max, emap_min)    # (0, 0, 1)
    i101 = energy_map_index([point1[0], point0[1], point1[2]], emap_max, emap_min)    # (1, 0, 1)
    i010 = energy_map_index([point0[0], point1[1], point0[2]], emap_max, emap_min)    # (0, 1, 0)
    i110 = energy_map_index([point1[0], point1[1], point0[2]], emap_max, emap_min)    # (1, 1, 0)
    i011 = energy_map_index([point0[0], point1[1], point1[2]], emap_max, emap_min)    # (0, 1, 1)
    i111 = energy_map_index(point1, emap_max, emap_min)                               # (1, 1, 1)

    c00 = emap[i000][atom_index] * (1 - dif[0]) + emap[i100][atom_index] * dif[0]
    c01 = emap[i001][atom_index] * (1 - dif[0]) + emap[i101][atom_index] * dif[0]
    c10 = emap[i010][atom_index] * (1 - dif[0]) + emap[i110][atom_index] * dif[0]
    c11 = emap[i011][atom_index] * (1 - dif[0]) + emap[i111][atom_index] * dif[0]

    c0 = c00 * (1 - dif[1]) + c10 * dif[1]
    c1 = c01 * (1 - dif[1]) + c11 * dif[1]

    c = c0 * (1 - dif[2]) + c1 * dif[2]

    return c


def run_interpenetration(sim_par, emap):
    """
    *** Not complete ***
    Run interpenetration algorithm with given simulation parameters and energy map.
    Returns simulation summary and structural information on the discovered structures.
    """
    # Initialize simulation parameters
    structure_energy_limit = sim_par['structure_energy_limit']
    atom_energy_limit = sim_par['atom_energy_limit']
    rotation_limit = sim_par['rotation_limit']
    rotation_freedom = sim_par['rotation_freedom']
    summary_percent = sim_par['summary_percent']

    # Initialize simulation variables
    Quat = Quaternion([0, 1, 1, 1])

    initial_coors = initial_coordinates(base_mof, emap, atom_list, 3E10)
    trial_limit = len(initial_coors) * rotation_limit
    rot_freedom = rot_freedom
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
                    x_angle = 2 * pi * math.floor(rand() * (rot_freedom)) / (rot_freedom)
                    y_angle = 2 * pi * math.floor(rand() * (rot_freedom)) / (rot_freedom)
                    z_angle = 2 * pi * math.floor(rand() * (rot_freedom)) / (rot_freedom)

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
                    structure['atom_coors'].append(first_point.xyz())  # Why first point not new_coor?
                    structure['pbc_coors'].append(new_coor.xyz())
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

                    point_energy = trilinear_interpolate(pbc_coor.xyz(), emap_atom_index, emap, emap_max, emap_min)
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


def check_extension(base_MOF, mobile_MOF, rotation_info, emap, emap_atom_list, energy_limit, ext_cut_off):
    """
    *** Not Complete ***
    Checks collision between interpenetrating layer and base layer for a determined distance.
    Distance is calculated from given ext_cut_off value which determines the packing amount of the
    interpenetrating layer.
    Each coordinate in the interpenetrating layer is checked for high energy values by applying
    perodic boundary conditions to the coordinate according to energy map of the base layer.
    """
    packing_factor = Packing.factor(mobile_MOF.uc_size, ext_cut_off)
    uc_vectors = Packing.uc_vectors(mobile_MOF.uc_size, mobile_MOF.uc_angle)
    trans_vec = Packing.translation_vectors(packing_factor, uc_vectors)
    packed_coors = Packing.uc_coors(trans_vec, packing_factor, uc_vectors, mobile_MOF.atom_coors)

    rotated_packed_coors = rotate_unit_cell(packed_coor, rotation_info)
    x_angle = rotation_info[0]
    y_angle = rotation_info[1]
    z_angle = rotation_info[2]

    collision = False
    for coor in rotated_packed_coors:

        if not collision:

            new_coor = Coor(coor)
            Q = Quaternion([1, new_coor.x, new_coor.y, new_coor.z])  # Might be a better way to do this
            Q = Quat.rotation(Q.xyz(), [0, 0, 0], [1, 0, 0], x_angle)
            Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 1, 0], y_angle)
            Q = Quat.rotation(Q.xyz(), [0, 0, 0], [0, 0, 1], z_angle)
            new_coor = Coor(Q.xyz())

            pbc_coor = new_coor.pbc(base_MOF.uc_size, base_MOF.uc_angle, base_MOF.frac_ucv)

            emap_index = energy_map_index(pbc_coor, emap_max, emap_min)
            atom_index = energy_map_atom_index(atom_name, emap_atom_list)

            energy = emap[emap_index][atom_index]
            if energy < energy_limit:
                continue
            else:
                collision = True
                break
        else:
            break

    return collision
