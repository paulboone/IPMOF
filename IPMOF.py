import os
import math

# Load 3rd party libraries
import yaml

# Load interpenetration python libraries
from ipmof.crystal import MOF
from ipmof.forcefield import read_ff_parameters
from ipmof.energymap import energy_map, get_mof_list, get_uniq_atom_list
from ipmof.interpenetration import run_interpenetration, check_extension, save_extension
from ipmof.core import core_mof_properties, core_mof_sort, core_mof_dir
# --------------------------------------------------------------------------------------------------
from ipmof.parameters import sim_dir_data as sim_dir
from ipmof.parameters import sim_par_data as sim_par

# Read excel file containing force field information
force_field = read_ff_parameters(sim_dir['force_field_path'], sim_par['force_field'])
# --------------------------------------------------------------------------------------------------
# Create MOf list from CoRE database
mof_properties = core_mof_properties(sim_dir['core_path'])

sorted_mofs = core_mof_sort(mof_properties, sort='void_fraction', limit=0.85)
mol2_list = core_mof_dir(sorted_mofs, sim_dir['mol2_dir'])

mof_list = get_mof_list(mol2_list, force_field)
print(mof_list)

for base_mof_index, base_mof in enumerate(mof_list):

    print('Base MOF selected as: ', base_mof.name)
    extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])

    print('Calculating emap for', base_mof.name, 'with atoms:', atom_list['atom'])

    # Calculate atom list for remaining MOFs
    atom_list = get_uniq_atom_list(mof_list[base_mof_index:])

    # Calculate energy map
    emap = energy_map(sim_par, base_mof, atom_list)

    for mobile_mof_index, mobile_mof in enumerate(mof_list):

        if mobile_mof_index >= base_mof_index:

            print('Mobile MOF selected as: ', mobile_mof.name)

            extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])

            # Run interpenetration
            summary, new_structures = run_interpenetration(sim_par, base_mof, mobile_mof, emap, atom_list)

            if len(new_structures) > 0:
                print('Interpenetration found. Structure count:', len(new_structures))
                print('MOF1:', base_mof.name, 'MOF2:', mobile_mof.name)
                # Get minimum energy structure by sorting total structure energies
                min_energy_structure = sorted(new_structures, key=lambda k: k['energy'])[0]

                # Add for loop here for min_energy_structures

                # Check for collision in the extended unitcell of new structure and energy map
                collision = check_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)
                print('Collision:', collision)

                # Get structure information for the interpenetrating structure
                ext_structure = save_extension(sim_par, base_mof, mobile_mof, emap, atom_list, min_energy_structure)

                # Extend MOF coordinates and get atom names and coordinates of extended unit cells of MOF object
                extended_structure = base_mof.extend_unit_cell(sim_par['cut_off'])

                # Create new MOF objects for base and mobile MOFs
                ext_base_mof = MOF(extended_structure, file_format='dict')
                ext_mobile_mof = MOF(ext_structure, file_format='dict')

                # Join base and mobile structure layers
                joined_mof = ext_base_mof.join(ext_mobile_mof, colorify=sim_par['export_colorify'])

                # Export to xyz format
                joined_mof.export(sim_dir['export_dir'], file_format=sim_par['export_format'])
            else:
                print('No interpenetration.')
