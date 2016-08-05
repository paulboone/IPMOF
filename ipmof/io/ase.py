# IPMOF file input/output library
# Date: August 2016
# Author: Kutay B. Sezginel
import os
from ase.io import read as ase_read
from ase.geometry import cell_to_cellpar as ase_cellpar


def read(file_path, input_format='cif'):
    # Read molecule file
    atoms = ase_read(file_path, format=input_format)

    # Get atom coordinates and names
    atom_coors = atoms.get_positions()
    atom_names = atoms.get_chemical_symbols()

    # Get unit cell parameters (converting from unit cell vectors)
    uc = ase_cellpar(atoms.cell)

    # Additional unit cell information (currently not used)
    sg = atoms.info['spacegroup'].symbol
    volume = atoms.get_volume()

    molecule = {'uc_size': [uc[0], uc[1], uc[2]],
                'uc_angle': [uc[3], uc[4], uc[5]],
                'atom_names': atom_names,
                'atom_coors': atom_coors}

    return molecule
