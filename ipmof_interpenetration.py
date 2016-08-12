import os
import sys

# Load IPMOF python libraries
from ipmof.crystal import MOF
from ipmof.forcefield import read_ff_parameters
from ipmof.energymap import import_energy_map, get_mof_list
from ipmof.interpenetration import run_interpenetration
from ipmof.core import core_mof_properties, core_mof_sort, core_mof_dir
from ipmof.parameters import export_init_txt
from ipmof.parameters import sim_dir_data as sim_dir    # Import simulation directories
from ipmof.parameters import sim_par_data as sim_par    # Import simulation parameters

# Read excel file containing force field information
force_field = read_ff_parameters(sim_dir['force_field_path'], sim_par['force_field'])

# Read MOF list from CoRE or from a given directory
if sim_par['core_database']:
    # Create MOf list from CoRE database
    mof_properties = core_mof_properties(sim_dir['core_path'])
    sorted_mofs = core_mof_sort(mof_properties, sort='void_fraction', limit=0.85)
    mol2_list = core_mof_dir(sorted_mofs, sim_dir['mof_dir'])
    mof_list = get_mof_list(mol2_list, force_field)
else:
    # Create MOF list by reading mol2 files from a directory
    mof_path_list = os.listdir(sim_dir['mof_dir'])
    mof_path_list = [os.path.join(sim_dir['mof_dir'], path) for path in mof_path_list]
    mof_list = get_mof_list(mof_path_list, force_field)

emap_name = os.listdir(sim_dir['energy_map_dir'])[0]
base_mof_dir = os.path.join(sim_dir['mof_dir'], emap_name.split('_emap')[0] + '.cif')
base_mof = MOF(base_mof_dir)
base_mof.force_field(force_field)
# Export initialization file containing MOF names and simulation parameters
export_init_txt(mof_list)

# Main Loop (Interpenetration)
print('-' * 80)
print('Energy map ->', base_mof.name)
print('-' * 80)
# Load energy map
atom_list, emap = import_energy_map(sim_par, sim_dir, base_mof.name)

for mobile_mof_index, mobile_mof in enumerate(mof_list):
    # Run interpenetration
    print('Running interpenetration for ', mobile_mof.name, '...')

    # Submit jobs here
    if sys.argv[-1] == 'q':
        # Load job server libraries
        from rq import Queue
        from redis import Redis
        import sjs

        # Load job queue
        sjs.load(os.path.join("settings", "sjs.yaml"))
        job_queue = sjs.get_job_queue()

        # Submit job
        job_queue.enqueue(run_interpenetration, base_mof, mobile_mof, emap, atom_list, sim_par, sim_dir)
    else:
        run_interpenetration(base_mof, mobile_mof, emap, atom_list, sim_par, sim_dir)
