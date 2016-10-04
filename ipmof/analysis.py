import os
import math
from ipmof.parameters import read_interpenetration_results, read_parameters
from ipmof.interpenetration import regenerate


def regenerate_structures(results_dir, colorify=True, file_format='cif', num=5):
    """
    Regenerates interenetrated structures from given results.yaml file.
    """
    results_path = os.path.join(results_dir, 'results.yaml')
    with open(results_path, 'r') as r:
        sim_par, structure_info, summary = read_interpenetration_results(results_path)

    if num == 'all':
        num = len(structure_info) - 1
    else:
        num = max(len(structure_info) - 1, num)

    s1_name = structure_info[0]['S1']
    s2_name = structure_info[0]['S2']
    export_dir = results_dir
    sim_par, sim_dir = read_parameters()

    for s_idx in range(1, num + 1):
        init_coor = structure_info[s_idx]['initial_coordinate']
        rot = structure_info[s_idx]['rotation']
        regenerate(s1_name, s2_name, rot, init_coor, sim_par, sim_dir,
                   export_dir, colorify=colorify, index=s_idx, format=file_format)
