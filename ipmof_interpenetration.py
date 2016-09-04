import os
import sys

# Load IPMOF python libraries
from ipmof.interpenetration import run_interpenetration, get_interpenetration_list
from ipmof.parameters import read_parameters

# Read simulation parameters and directories
sim_par, sim_dir = read_parameters()

# Get list of interpenetrating MOFs
interpenetration_list = get_interpenetration_list(sim_par, sim_dir)
print('Initializing interpenetration for', len(interpenetration_list), 'MOF combinations...')

for ip_index, interpenetration_path in enumerate(interpenetration_list, start=1):

    emap_path, emap_mof_path, ip_mof_path = interpenetration_path
    print('-' * 80 + '\n' + str(ip_index), 'Energy map ->', os.path.basename(emap_path))
    print(' ' * len(str(ip_index)), 'Interpenetration ->', os.path.basename(ip_mof_path) + '\n' + '-' * 80)

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
        job_queue.enqueue(run_interpenetration, interpenetration_path, sim_par, sim_dir)
    else:
        run_interpenetration(interpenetration_path, sim_par, sim_dir)
