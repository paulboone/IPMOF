import os
import math

# Load 3rd party libraries
import yaml

# Load interpenetration python libraries
from ipmof.crystal import MOF
from ipmof.forcefield import read_ff_parameters
from ipmof.energymap import energy_map, get_mof_list, get_mof_file_list, get_uniq_atom_list
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

# Export initialization file containing MOF names and simulation parameters
export_init_txt(mof_list)
# --------------------------------------------------------------------------------------------------
# Main Loop (Energy Map) ---------------------------------------------------------------------------
for base_mof_index, base_mof in enumerate(mof_list):

    print('-----------------------------------------------------------------------------------')
    print(base_mof_index, 'Base ->', base_mof.name)
    print('-----------------------------------------------------------------------------------')

    # Calculate atom list for remaining MOFs
    atom_list = get_uniq_atom_list(mof_list[base_mof_index:])

    # Extend base MOF for energy map calculation
    extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])

    # Calculate energy map -------------------------------------------------------------------------
    emap = energy_map(sim_par, base_mof, atom_list)

    for mobile_mof_index, mobile_mof in enumerate(mof_list):

        if mobile_mof_index >= base_mof_index:

            # extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])

            # Run interpenetration -----------------------------------------------------------------
            summary, new_structures = run_interpenetration(sim_par, base_mof, mobile_mof, emap, atom_list)

            # Create export directory and export summary -------------------------------------------
            export_dir = os.path.join(sim_dir['export_dir'], base_mof.name + '_' + mobile_mof.name)
            if not os.path.isdir(export_dir):
                os.mkdir(export_dir)
            if sim_par['export_summary']:
                export_summary_txt(export_dir, summary, base_mof, mobile_mof)

            # Export Min Energy Structures ---------------------------------------------------------
            if len(new_structures) > 0:
                print(mobile_mof_index, base_mof.name, '--', mobile_mof.name, '-> (+) Structure:', len(new_structures))

                export_count = min(len(new_structures), sim_par['export_structures'])
                for export_index in range(export_count):
                    # Get minimum energy structure by sorting total structure energies
                    min_energy_structure = sorted(new_structures, key=lambda k: k['energy'])[export_index]

                    # Check for collision in the extended unitcell of new structure and energy map
                    collision = check_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)

                    # Print structure information --------------------------------------------------
                    structure_info = '+ Structure ' + str(export_index) + ' -> Collision: '
                    structure_info = str(collision) + ' | Energy: '
                    structure_info += str(round(min_energy_structure['energy'], 4))
                    structure_info += '\n              -> Rotation: '
                    structure_info += str(min_energy_structure['rotation']) + ' First Point: '
                    structure_info += str(min_energy_structure['first_point']) + '\n'
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
                        packed_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
                        packed_mobile_mof = MOF(packed_structure, file_format='dict')
                        packed_base_mof = MOF(extended_structure, file_format='dict')
                        joined_packed_mof = packed_base_mof.join(packed_mobile_mof, colorify=False)
                        joined_packed_mof.export(export_dir, file_format=sim_par['export_format'])

                    if sim_par['export_packed_color']:
                        # Pack new structure by using rotation and first point information
                        packed_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
                        packed_mobile_mof = MOF(packed_structure, file_format='dict')
                        packed_base_mof = MOF(extended_structure, file_format='dict')
                        joined_packed_mof = packed_base_mof.join(packed_mobile_mof, colorify=True)
                        joined_packed_mof.export(export_dir, file_format=sim_par['export_format'])
            else:
                print(mobile_mof_index, base_mof.name, '--', mobile_mof.name, '-> (-)')
