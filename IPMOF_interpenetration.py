import os
import math

# Load 3rd party libraries
import yaml
import numpy as np

# Load interpenetration python libraries
from ipmof.crystal import MOF
from ipmof.forcefield import read_ff_parameters
from ipmof.energymap import import_energy_map, get_mof_list, get_mof_file_list, get_uniq_atom_list
from ipmof.interpenetration import run_interpenetration, check_extension, save_extension
from ipmof.core import core_mof_properties, core_mof_sort, core_mof_dir
from ipmof.parameters import export_init_txt, export_summary_txt
from ipmof.parameters import sim_dir_data as sim_dir    # Import simulation directories
from ipmof.parameters import sim_par_data as sim_par    # Import simulation parameters
# --------------------------------------------------------------------------------------------------
# Read excel file containing force field information
force_field = read_ff_parameters(sim_dir['force_field_path'], sim_par['force_field'])

# Read MOF list from CoRE or from a given directory
if sim_par['core_database']:
    # Create MOf list from CoRE database
    mof_properties = core_mof_properties(sim_dir['core_path'])
    sorted_mofs = core_mof_sort(mof_properties, sort='void_fraction', limit=0.85)
    mol2_list = core_mof_dir(sorted_mofs, sim_dir['mol2_dir'])
    mof_list = get_mof_list(mol2_list, force_field)
else:
    # Create MOF list by reading mol2 files from a directory
    mof_list = get_mof_file_list(sim_dir['mol2_dir'], 'mol2', force_field)

emap_name = os.listdir(sim_dir['energy_map_dir'])[0]
base_mof_dir = os.path.join(sim_dir['mol2_dir'], emap_name.split('_')[0] + '.mol2')
base_mof = MOF(base_mof_dir)
base_mof.force_field(force_field)

# Export initialization file containing MOF names and simulation parameters
export_init_txt(mof_list)
# --------------------------------------------------------------------------------------------------
# Main Loop (Interpenetration) ---------------------------------------------------------------------
print('-----------------------------------------------------------------------------------')
print('Energy map ->', base_mof.name)
print('-----------------------------------------------------------------------------------')
# Load energy map ------------------------------------------------------------------------------
atom_list, emap = import_energy_map(sim_par, sim_dir, base_mof.name)

for mobile_mof_index, mobile_mof in enumerate(mof_list):
    # Run interpenetration -----------------------------------------------------------------
    print('Running interpenetration for ', mobile_mof.name, '...')
    summary, new_structures = run_interpenetration(sim_par, base_mof, mobile_mof, emap, atom_list)

    # Create export directory and export summary -------------------------------------------
    export_dir = os.path.join(sim_dir['export_dir'], base_mof.name + '_' + mobile_mof.name)
    if not os.path.isdir(export_dir):
        os.mkdir(export_dir)
    if sim_par['export_summary']:
        export_summary_txt(export_dir, summary, base_mof, mobile_mof)

    # Export Min Energy Structures ---------------------------------------------------------
    if len(new_structures) > 0:
        print(base_mof.name, '--', mobile_mof.name, '-> (+) Structure:', len(new_structures))

        export_count = min(len(new_structures), sim_par['export_structures'])
        for export_index in range(export_count):
            # Get minimum energy structure by sorting total structure energies
            min_energy_structure = sorted(new_structures, key=lambda k: k['energy'])[export_index]

            # Check for collision in the extended unitcell of new structure and energy map
            collision = check_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)

            # Print structure information --------------------------------------------------
            rot_x = str(round(math.degrees(min_energy_structure['rotation'][0])))
            rot_y = str(round(math.degrees(min_energy_structure['rotation'][1])))
            rot_z = str(round(math.degrees(min_energy_structure['rotation'][2])))
            fp_x = str(round(min_energy_structure['first_point'][0], 3))
            fp_y = str(round(min_energy_structure['first_point'][1], 3))
            fp_z = str(round(min_energy_structure['first_point'][2], 3))
            structure_info = '\tEnergy: ' + str(round(min_energy_structure['energy'], 2))
            structure_info += ' | Collision: ' + str(collision) + '\n'
            structure_info += '\tRotation x: ' + rot_x + ' y: ' + rot_y + ' z: ' + rot_z
            structure_info += ' | First Point x: ' + fp_x + ' y: ' + fp_y + ' z: ' + fp_z
            print(structure_info)

            # Record new structure ---------------------------------------------------------
            new_structure = {'atom_names': min_energy_structure['atom_names'], 'name': mobile_mof.name}
            if sim_par['export_pbc']:
                new_structure['atom_coors'] = min_energy_structure['pbc_coors']
            else:
                new_structure['atom_coors'] = min_energy_structure['atom_coors']
            new_mobile_mof = MOF(new_structure, file_format='dict')

            # Export structures ------------------------------------------------------------
            if sim_par['export_single']:
                new_mobile_mof = MOF(new_structure, file_format='dict')
                joined_mof = base_mof.join(new_mobile_mof, colorify=False)
                joined_mof.name += str(export_index)
                joined_mof.export(export_dir, file_format=sim_par['export_format'])

            if sim_par['export_single_color']:
                # Join base and mobile structure layers
                new_mobile_mof = MOF(new_structure, file_format='dict')
                joined_mof_color = base_mof.join(new_mobile_mof, colorify=True)
                joined_mof_color.name += str(export_index) + 'C'
                joined_mof_color.export(export_dir, file_format=sim_par['export_format'])

            if sim_par['export_packed']:
                # Pack new structure by using rotation and first point information
                extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])
                packed_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
                packed_mobile_mof = MOF(packed_structure, file_format='dict')
                packed_base_mof = MOF(extended_structure, file_format='dict')
                joined_packed_mof = packed_base_mof.join(packed_mobile_mof, colorify=False)
                joined_packed_mof.export(export_dir, file_format=sim_par['export_format'])

            if sim_par['export_packed_color']:
                # Pack new structure by using rotation and first point information
                extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])
                packed_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
                packed_mobile_mof = MOF(packed_structure, file_format='dict')
                packed_base_mof = MOF(extended_structure, file_format='dict')
                joined_packed_mof = packed_base_mof.join(packed_mobile_mof, colorify=True)
                joined_packed_mof.export(export_dir, file_format=sim_par['export_format'])
    else:
        print(mobile_mof_index, base_mof.name, '--', mobile_mof.name, '-> (-)')
