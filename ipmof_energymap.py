import os
import sys

# Load IPMOF python libraries
from ipmof.forcefield import read_ff_parameters
from ipmof.energymap import energy_map, get_mof_list, get_mof_file_list, get_uniq_atom_list
from ipmof.core import core_mof_properties, core_mof_sort, core_mof_dir
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

# Calculate atom list for remaining MOFs
if sim_par['energy_map_atom_list'] == 'uniq':
    atom_list = get_uniq_atom_list(mof_list)
elif sim_par['energy_map_atom_list'] == 'full':
    atom_list = force_field

# Export initialization file containing MOF names and simulation parameters
print('Starting energy map calculation with grid size:', sim_par['grid_size'],
      'and cut off radius:', sim_par['cut_off'])
print('Atom list ->', atom_list['atom'])
print('Energy map(s) will be exported in', sim_par['energy_map_type'], 'format')
# --------------------------------------------------------------------------------------------------
# Main Loop (Energy Map) ---------------------------------------------------------------------------
for mof_index, mof in enumerate(mof_list):

    # Extend base MOF for energy map calculation
    extended_structure = mof.extend_unit_cell(sim_par['cut_off'])

    print('-----------------------------------------------------------------------------------')
    print(mof_index, 'Calculating energy map for ->', mof.name)

    # Submit jobs here
    if sys.argv[-1] == 'q':
        # Load job server libraries
        from rq import Queue
        from redis import Redis
        import sjs

        # Load job queue
        sjs.load(os.path.join("settings", "sjs.yaml"))
        job_queue = sjs.get_job_queue()

        # Calculate energy map ---------------------------------------------------------------------
        job_queue.enqueue(energy_map, sim_par, mof, atom_list)
    else:
        emap = energy_map(sim_par, mof, atom_list)
