# IPMOF file input/output library
# Date: August 2016
# Author: Kutay B. Sezginel
import os
import openbabel as ob


def read(file_path, input_format='cif'):
    """
    Read file from a given path.
    """
    conv = ob.OBConversion()
    conv.SetInFormat(input_format)
    mol = ob.OBMol()
    conv.ReadFile(mol, file_path)

    # Get atom names and coordinates
    atom_names = []
    atom_coors = []
    for atom_index in range(1, mol.NumAtoms() + 1):
        atom = mol.GetAtom(atom_index)
        atom_name = str(atom.GetType())
        # Remove digits from atom name
        atom_name = ''.join([char for char in atom_name if not char.isdigit()])
        atom_names.append(atom_name)
        atom_coors.append([atom.x(), atom.y(), atom.z()])

    # Get unit cell information
    uc = ob.toUnitCell(mol.GetData(ob.UnitCell))

    # Additional unit cell information (currently not used)
    sg = uc.GetSpaceGroupName()
    volume = uc.GetCellVolume()

    molecule = {'uc_size': [uc.GetA(), uc.GetB(), uc.GetC()],
                'uc_angle': [uc.GetAlpha(), uc.GetBeta(), uc.GetGamma()],
                'atom_names': atom_names,
                'atom_coors': atom_coors}

    return molecule


def write(mol, export_path, output_format='xyz'):
    conv = ob.OBConversion()
    conv.SetOutFormat(output_format)
    conv.WriteFile(mol, export_path + '.' + output_format)
