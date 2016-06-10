import os
import math
import yaml
import sys
# --------------------------------------------------------------------------------------------------
input_dir = os.getcwd()

sim_par_path = os.path.join(input_dir, 'sim_par.yaml')
sim_dir_path = os.path.join(input_dir, 'sim_dir_linux.yaml')

# Read sim par yaml file
sim_par = yaml.load(open(sim_par_path, 'r'))
sim_dir = yaml.load(open(sim_dir_path, 'r'))

# Load interpenetration python libraries
sys.path.append(sim_dir['python_lib_dir'])
from forcefield import read_ff_parameters
from crystal import MOF, extend_unit_cell, export_xyz
from energymap import energy_map, get_mof_list, get_uniq_atom_list
from interpenetration import run_interpenetration, check_extension, save_extension, join_structures

# Read excel file containing force field information
force_field = read_ff_parameters(sim_dir['excel_file_path'], sim_par['force_field'])
# --------------------------------------------------------------------------------------------------
# Create list of MOFs
mol2_list = get_mof_list(sim_dir['mol2_dir'], '.mol2')
print(mol2_list)

mof_list = [mol2_list[4]]

for base_index, base_mof_selection in enumerate(mof_list):
    for mobile_index, mobile_mof_selection in enumerate(mof_list):

        if mobile_index >= base_index:
            # Read mol2 files and initialize MOF objects
            mol2_path = os.path.join(sim_dir['mol2_dir'], mol2_list[base_mof_index])
            base_mof = MOF(mol2_path)
            base_mof.force_field(force_field)
            print('Base MOF selected as: ', base_mof.name)

            mol2_path = os.path.join(sim_dir['mol2_dir'], mol2_list[mobile_mof_index])
            mobile_mof = MOF(mol2_path)
            mobile_mof.force_field(force_field)
            print('Mobile MOF selected as: ', mobile_mof.name)

            extended_structure = extend_unit_cell(base_mof, sim_par['cut_off'])

            atom_list = get_uniq_atom_list([mobile_mof])
            print('Calculating emap for', base_mof.name, 'with atoms:', atom_list['atom'])

            # Calculate energy map
            emap = energy_map(sim_par, base_mof, atom_list)

            # Run interpenetration
            summary, new_structures = run_interpenetration(sim_par, base_mof, mobile_mof, emap, atom_list)

            if len(new_structures) > 0:
                print('Interpenetration found. Structure count:', len(new_structures))
                # Get minimum energy structure by sorting total structure energies
                min_energy_structure = sorted(new_structures, key=lambda k: k['energy'])[0]

                # Check for collision in the extended unitcell of new structure and energy map
                collision = check_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
                print('Collision:', collision)

                # Get structure information for the interpenetrating structure
                ext_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)

                # Extend MOF coordinates and get atom names and coordinates of extended unit cells of MOF object
                extended_structure = extend_unit_cell(base_mof, sim_par['ext_cut_off'])

                # Join base and mobile structure layers
                joined_structure = join_structures(extended_structure, ext_structure, colorify=True)

                # Export to xyz format
                structure_name = base_mof.name + '_' + mobile_mof.name
                export_xyz(joined_structure['atom_coors'], joined_structure['atom_names'], structure_name, sim_dir['export_dir'])
            else:
                print('No interpenetration.')
