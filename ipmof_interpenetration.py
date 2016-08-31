import os
import sys

# Load IPMOF python libraries
from ipmof.interpenetration import run_interpenetration, get_interpenetration_list
from ipmof.parameters import read_parameters

# Read simulation parameters and directories
sim_par, sim_dir = read_parameters()

# Get list of interpenetrating MOFs
interpenetration_list = get_interpenetration_list(sim_par, sim_dir)
print('Initializing interpenetration for', len(interpenetration_list['emap_path']), 'MOF combinations...')

for list_index, emap_path in enumerate(interpenetration_list['emap_path']):
    ip_mof_path = interpenetration_list['ip_mof_path'][list_index]
    emap_mof_path = interpenetration_list['emap_mof_path'][list_index]

    ip_index = str(list_index + 1)
    print('-' * 80 + '\n' + ip_index, 'Energy map ->', os.path.basename(emap_path))
    print(' ' * len(ip_index), 'Interpenetration ->', os.path.basename(ip_mof_path) + '\n' + '-' * 80)

    # Run interpenetration
    if sys.argv[-1] == 'q':
        # Load job server libraries
        from rq import Queue
        from redis import Redis
        import sjs

        # Load job queue
        sjs.load(os.path.join("settings", "sjs.yaml"))
        job_queue = sjs.get_job_queue()

        # Run interpenetration
        job_queue.enqueue(run_interpenetration, emap_mof_path, ip_mof_path, emap_path, sim_par, sim_dir)
    else:
        run_interpenetration(emap_mof_path, ip_mof_path, emap_path, sim_par, sim_dir)
