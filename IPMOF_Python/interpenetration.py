# Quaternion class
# - Multiplication  - Division  - Inverse   - Rotation  - Output List <Q.xyz()>
#  >>> Q = Quaternion([0, 1, 1, 1])
# Date: June 2016
# Author: Kutay B. Sezginel
import math
from crystal import Packing, MOF
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


def check_extension(MOF, emap, ext_cut_off):
    """
    *** Not Complete ***
    Checks collision between interpenetrating layer and base layer for a determined distance.
    Distance is calculated from given ext_cut_off value which determines the packing amount of the
    interpenetrating layer.
    Each coordinate in the interpenetrating layer is checked for high energy values by applying
    perodic boundary conditions to the coordinate according to energy map of the base layer.
    """
    packing_factor = Packing.factor(MOF.uc_size, ext_cut_off)
    uc_vectors = Packing.uc_vectors(MOF.uc_size, MOF.uc_angle)
    trans_vec = Packing.translation_vectors(packing_factor, uc_vectors)
    MOF.packed_coors = Packing.uc_coors(trans_vec, packing_factor, uc_vectors, MOF.atom_coors)

    rotated_packed_coors = rotate_unit_cell(MOF)

    for coor in rotated_packed_coors:
        atom_coor = Coor(coor)
        pbc_coor = pbc_coor.pbc(uc_size, uc_angle)

        emap_index = energy_map_index(pbc_coor)
        atom_index = find_atom_index(emap)

        energy = emap[emap_index][atom_index]
        if energy < energy_limit:
            continue
        else:
            break