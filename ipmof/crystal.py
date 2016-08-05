# IPMOF Crystal Functions
# Date: June 2016
# Author: Kutay B. Sezginel
import math
import os
from ipmof.forcefield import get_ff_parameters
import ipmof.io as io


class MOF:
    """
    MOF class that holds coordinate, atom name, and unit cell information
    """
    def __init__(self, file_path, file_format='cif', reader='ase'):
        """
        Initialize MOF name, unit cell volume and parameters, atom names and coordinates, and
        unique atom names and coordinates.
        reader -> ('babel' / 'ase')
        file_format -> ('cif' / 'dict' / 'mol2' / ...)
        """
        if file_format == 'dict':
            molecule = file_path
            self.atom_coors = molecule['atom_coors']
            self.atom_names = molecule['atom_names']
            self.name = molecule['name']
        else:
            self.path = file_path
            self.name, fileformat = os.path.basename(file_path).split('.')
            # Import reader library from ipmof.io and read structure file into a dictionary
            reader_lib = __import__('ipmof.io.' + reader, fromlist=[''])
            molecule = reader_lib.read(file_path, input_format=file_format)

            self.atom_coors = molecule['atom_coors']
            self.atom_names = molecule['atom_names']
            self.uc_size = molecule['uc_size']
            self.uc_angle = molecule['uc_angle']
            self.separate_atoms()
            self.unit_cell_volume()

    def __repr__(self):
        return "<MOF object: %s>" % (self.name)

    def __str__(self):
        return self.name

    def __len__(self):
        return len(self.atom_coors)

    def force_field(self, ff_type):
        """
        Initializes force field parameters according to unique atom names.
        """
        ff_param = get_ff_parameters(self.uniq_atom_names, ff_type)
        self.sigma = []
        self.epsilon = []
        for i in range(len(ff_param)):
            self.uniq_atom_names[i] = ff_param[i][0]
            self.sigma.append(ff_param[i][1])
            self.epsilon.append(ff_param[i][2])

    def extend_unit_cell(self, cut_off):
        """
        Extends unit cell of a MOF object according to a given cut_off value.
        The structure is extended so that all the points in the center unit cell are
        at least *cut_off Angstroms away.
        This enables energy map calculation by providing coordinates for surrounding atoms.
        Also enables checking for collisions in the extended unit cell of mobile layer
        and energy map.
        The *ext_cut_off parameter in the sim_par input file determines the amount of packing.
        - A high cut off value such as 100 Angstrom can be used to ensure there is no collisions
        between the interpenetrating layers.
        """
        self.packing_factor = Packing.factor(self.uc_size, cut_off)
        uc_vectors = Packing.uc_vectors(self.uc_size, self.uc_angle)
        trans_vec = Packing.translation_vectors(self.packing_factor, uc_vectors)
        self.packed_coors = Packing.uc_coors(trans_vec, self.packing_factor, uc_vectors, self.atom_coors)
        self.edge_points = Packing.edge_points(uc_vectors)

        extended_structure = {'atom_names': [], 'atom_coors': [], 'name': self.name}
        for unit_cell in self.packed_coors:
            for coor_index, coor in enumerate(unit_cell):
                atom_name = self.atom_names[coor_index]
                extended_structure['atom_names'].append(atom_name)
                extended_structure['atom_coors'].append(coor)

        return extended_structure

    def separate_atoms(self):
        """
        Identifies different types of atoms in given atom coordinates and atom names
        Returns unique atom coordinates and unique atom names
        """
        uniq_atom_names = []
        for atom in self.atom_names:
            if atom not in uniq_atom_names:
                uniq_atom_names.append(atom)

        uniq_atom_coors = [[] for i in range(len(uniq_atom_names))]
        for atom_index in range(len(self.atom_coors)):
            for uniq_atom_index, uniq_atom in enumerate(uniq_atom_names):
                if self.atom_names[atom_index] == uniq_atom:
                    uniq_atom_coors[uniq_atom_index].append(self.atom_coors[atom_index])

        self.uniq_atom_names = uniq_atom_names
        self.uniq_atom_coors = uniq_atom_coors

    def calculate_cut_off(self):
        """
        Calculate cut-off radius as Rc = L/2 from a given MOF object.
        """
        width_a = self.ucv / (self.uc_size[1] * self.uc_size[2] / math.sin(math.radians(self.uc_angle[0])))
        width_b = self.ucv / (self.uc_size[0] * self.uc_size[2] / math.sin(math.radians(self.uc_angle[1])))
        width_c = self.ucv / (self.uc_size[0] * self.uc_size[1] / math.sin(math.radians(self.uc_angle[2])))
        self.cut_off = min(width_a / 2, width_b / 2, width_c / 2)

    def unit_cell_volume(self):
        """
        Calculates unit cell volume of a given MOF object.
        """
        a = self.uc_size[0]
        b = self.uc_size[1]
        c = self.uc_size[2]
        alp = math.radians(self.uc_angle[0])
        bet = math.radians(self.uc_angle[1])
        gam = math.radians(self.uc_angle[2])

        volume = 1 - math.cos(alp)**2 - math.cos(bet)**2 - math.cos(gam)**2
        volume += 2 * math.cos(alp) * math.cos(bet) * math.cos(gam)
        volume = a * b * c * math.sqrt(volume)
        frac_volume = volume / (a * b * c)

        self.ucv = volume
        self.frac_ucv = frac_volume

    def join(self, new_mof, colorify=False, atom_color=['C', 'O']):
        """
        Combines atom names and coordinates of two given structure dictionaries.
        The structure dictionaries should be in following format:
         >>> structure = {'atom_names': [*atom names], 'atom_coors':[*atom coors]}
        """
        base_atom_name = atom_color[0]
        new_atom_name = atom_color[1]
        joined_atom_names = []
        joined_atom_coors = []

        if colorify:
            for coor in new_mof.atom_coors:
                joined_atom_names.append(new_atom_name)
                joined_atom_coors.append(coor)

            for coor in self.atom_coors:
                joined_atom_names.append(base_atom_name)
                joined_atom_coors.append(coor)

        else:
            for atom, coor in zip(new_mof.atom_names, new_mof.atom_coors):
                joined_atom_names.append(atom)
                joined_atom_coors.append(coor)

            for atom, coor in zip(self.atom_names, self.atom_coors):
                joined_atom_names.append(atom)
                joined_atom_coors.append(coor)

        joined_structure = {'atom_names': joined_atom_names, 'atom_coors': joined_atom_coors}
        joined_structure['name'] = self.name + '_' + new_mof.name
        joined_mof = MOF(joined_structure, file_format='dict')

        return joined_mof

    def export(self, export_dir, file_format='xyz'):
        """
        Export MOF atom coordinates and names in .xyz format.

        Example usage:
         >>> mof.export(export_dir, file_format='xyz')
        """
        if file_format == 'xyz':
            xyz_path = os.path.join(export_dir, self.name + '.xyz')
            if os.path.exists(xyz_path):
                os.remove(xyz_path)
            xyz_file = open(xyz_path, 'w')
            xyz_file.write(str(len(self.atom_coors)) + '\n')
            xyz_file.write(self.name + '\n')

            for atom, coor in zip(self.atom_names, self.atom_coors):
                xyz_file.write(atom + ' ' + str(coor[0]) + ' ' + str(coor[1]) + ' ' + str(coor[2]) + '\n')
            xyz_file.close()


class Packing:
    """
    Packing class containing functions used for packing unit cells
    """
    @classmethod
    def factor(cls, uc_size, cut_off):
        """
        Calculate packing factor for given unit cell size and cut off radius
        """
        packing_factor = []
        for uc_length in uc_size:
            packing_factor.append(math.ceil(cut_off / uc_length) * 2 + 1)
        return packing_factor

    @classmethod
    def uc_vectors(cls, uc_size, uc_angle):
        """
        Calculate unit cell vectors for given unit cell size and angles
        """
        a = uc_size[0]
        b = uc_size[1]
        c = uc_size[2]
        alpha = math.radians(uc_angle[0])
        beta = math.radians(uc_angle[1])
        gamma = math.radians(uc_angle[2])

        x_v = [a, 0, 0]
        y_v = [b * math.cos(gamma), b * math.sin(gamma), 0]
        z_v = [0.0] * 3
        z_v[0] = c * math.cos(beta)
        z_v[1] = (c * b * math.cos(alpha) - y_v[0] * z_v[0]) / y_v[1]
        z_v[2] = math.sqrt(c * c - z_v[0] * z_v[0] - z_v[1] * z_v[1])
        uc_vectors = [x_v, y_v, z_v]
        return uc_vectors

    @classmethod
    def translation_vectors(cls, packing_factor, uc_vectors):
        """
        Calculate translation vectors for given packing factor and uc vectors
        """
        packing_amount = []
        for x in range(packing_factor[0]):
            for y in range(packing_factor[1]):
                for z in range(packing_factor[2]):
                    packing_amount.append([x, y, z])

        x_v = uc_vectors[0]
        y_v = uc_vectors[1]
        z_v = uc_vectors[2]
        translation_vectors = []
        for pack in packing_amount:
            x_trans = x_v[0] * pack[0] + y_v[0] * pack[1] + z_v[0] * pack[2]
            y_trans = x_v[1] * pack[0] + y_v[1] * pack[1] + z_v[1] * pack[2]
            z_trans = x_v[2] * pack[0] + y_v[2] * pack[1] + z_v[2] * pack[2]
            translation_vectors.append([x_trans, y_trans, z_trans])

        return translation_vectors

    @classmethod
    def uc_coors(cls, translation_vectors, packing_factors, uc_vectors, atom_coors):
        """
        Calculate packed coordinates for given:
        - translation vectors  - packing factor     - unit cell vectors    - atom coordinates
        """
        x_v = uc_vectors[0]
        y_v = uc_vectors[1]
        z_v = uc_vectors[2]

        packed_coors = [[] for i in range(len(translation_vectors))]
        translation_factor = []
        origin_trans_vec = []
        for dim, factor in enumerate(packing_factors):
            translation_factor.append((factor - 1) / 2)
            origin_translation = (factor - 1) / 2 * x_v[dim] + (factor - 1) / 2 * y_v[dim]
            origin_translation += (factor - 1) / 2 * z_v[dim]
            origin_trans_vec.append(origin_translation)

        packing_index = 0
        for translation in translation_vectors:
            for coor in atom_coors:
                x = coor[0] + translation[0] - origin_trans_vec[0]
                y = coor[1] + translation[1] - origin_trans_vec[1]
                z = coor[2] + translation[2] - origin_trans_vec[2]
                packed_coors[packing_index].append([x, y, z])
            packing_index += 1

        return packed_coors

    @classmethod
    def edge_points(cls, uc_vectors):
        """
        Calculate coordinates of unit cell edges in order of:
        (0, 0, 0) - (a, 0, 0) - (b, 0, 0) - (c, 0, 0)
        (a, b, 0) - (0, b, c) - (a, 0, c) - (a, b, c)
        """
        uc_edges = []
        uc_edges.append([0, 0, 0])

        # (a, 0, 0) - (b, 0, 0) - (c, 0, 0)
        for vec in uc_vectors:
            uc_edges.append(vec)

        # (a, b, 0)
        uc_edges.append([uc_vectors[0][0] + uc_vectors[1][0], uc_vectors[0][1] + uc_vectors[1][1],
                        uc_vectors[0][2] + uc_vectors[1][2]])
        # (0, b, c)
        uc_edges.append([uc_vectors[1][0] + uc_vectors[2][0], uc_vectors[1][1] + uc_vectors[2][1],
                        uc_vectors[1][2] + uc_vectors[2][2]])
        # (a, 0, c)
        uc_edges.append([uc_vectors[0][0] + uc_vectors[2][0], uc_vectors[0][1] + uc_vectors[2][1],
                        uc_vectors[0][2] + uc_vectors[2][2]])
        # (a, b, c)
        uc_edges.append([uc_vectors[0][0] + uc_vectors[1][0] + uc_vectors[2][0],
                        uc_vectors[0][1] + uc_vectors[1][1] + uc_vectors[2][1],
                        uc_vectors[0][2] + uc_vectors[1][2] + uc_vectors[2][2]])
        return uc_edges
