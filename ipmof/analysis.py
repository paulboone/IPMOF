import os
import math
from ipmof.parameters import read_interpenetration_results, read_parameters
from ipmof.interpenetration import regenerate
from tabulate import tabulate


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
        num = min(len(structure_info) - 1, num)

    s1_name = structure_info[0]['S1']
    s2_name = structure_info[0]['S2']
    export_dir = results_dir
    sim_par, sim_dir = read_parameters()

    for s_idx in range(1, num + 1):
        init_coor = structure_info[s_idx]['initial_coordinate']
        rot = structure_info[s_idx]['rotation']
        regenerate(s1_name, s2_name, rot, init_coor, sim_par, sim_dir,
                   export_dir, colorify=colorify, index=s_idx, format=file_format)


def summarize_results(results_dir, table=True, full=False, sortby='structure'):
    """
    Summarize interpenetration results in a given directory and export summary file.
    Optional arguments:
     - table (True or False):
        generate table for min. energy of each MOF combination

     - full (True or False):
        generate table for all combinations (True) or only for combinations
        that resulted in interpenetration (False)

     - sortby ('structure', 'energy', 'collision'):
        sort the table according to number of structures discovered,
        energy of the structure or extended collision
    """
    total_num_of_structures = 0
    table_lines = []
    no_results = []
    ip_count = dict(homo=0, hetero=0, homo_total=0, hetero_total=0, no_res=0, no_structure=0)
    homo_structures = 0
    hetero_structures = 0

    for mof_combination in os.listdir(results_dir):
        results_path = os.path.join(results_dir, mof_combination, 'results.yaml')

        if not os.path.isfile(results_path):
            no_results.append(mof_combination)
            ip_count['no_res'] += 1
        else:
            sim_par, structure_info, summary = read_interpenetration_results(results_path)
            mof_names = structure_info[0]['S1'] + ' + ' + structure_info[0]['S2']
            num_of_structures = structure_info[0]['Structures']
            total_num_of_structures += num_of_structures

            if structure_info[0]['S1'] == structure_info[0]['S2']:
                ip_count['homo_total'] += 1
            else:
                ip_count['hetero_total'] += 1

            if num_of_structures > 0:
                min_energy = structure_info[1]['energy']
                collision = structure_info[1]['collision']['exist']
                rotation = structure_info[1]['rotation']
                if structure_info[0]['S1'] == structure_info[0]['S2']:
                    ip_count['homo'] += 1
                    homo_structures += num_of_structures
                else:
                    ip_count['hetero'] += 1
                    hetero_structures += num_of_structures
                if table and not full:
                    table_lines.append([mof_names, num_of_structures, min_energy, collision, rotation])
            else:
                ip_count['no_structure'] += 1
                min_energy = collision = rotation = ' - - - '
            if table and full:
                table_lines.append([mof_names, num_of_structures, min_energy, collision, rotation])

    ip_count['total'] = ip_count['hetero_total'] + ip_count['homo_total']
    ip_count['structure'] = ip_count['hetero'] + ip_count['homo']

    sdr = round(ip_count['structure'] / ip_count['total'] * 100, 2)
    sdr_homo_t = round(ip_count['homo'] / ip_count['total'] * 100, 2)
    sdr_hetero_t = round(ip_count['hetero'] / ip_count['total'] * 100, 2)
    sdpt = round(total_num_of_structures / ip_count['total'], 2)

    if ip_count['homo_total'] == 0:
        sdr_homo_i = '---'
        sdpt_homo_t = '---'
    else:
        sdr_homo_i = round(ip_count['homo'] / ip_count['homo_total'] * 100, 2)
        sdpt_homo_t = round(homo_structures / ip_count['homo_total'], 2)

    if ip_count['hetero_total'] == 0:
        sdr_hetero_i = '---'
        sdpt_hetero_t = '---'
    else:
        sdr_hetero_i = round(ip_count['hetero'] / ip_count['hetero_total'] * 100, 2)
        sdpt_hetero_t = round(hetero_structures / ip_count['hetero_total'], 2)

    if ip_count['hetero'] == 0:
        sdpt_hetero_i = '---'
    else:
        sdpt_hetero_i = round(hetero_structures / ip_count['hetero'], 2)

    if ip_count['homo'] == 0:
        sdpt_homo_i = '---'
    else:
        sdpt_homo_i = round(homo_structures / ip_count['homo'], 2)

    analysis_text = 'Total MOF combinations:  ' + str(ip_count['total'])
    analysis_text += '\tFailed jobs:  ' + str(ip_count['no_res']) + '\n'
    analysis_text += '\tHetero:  ' + str(ip_count['hetero_total'])
    analysis_text += '\tHomo:  ' + str(ip_count['homo_total']) + '\n'
    analysis_text += '\nMOF combinations with structures:  ' + str(ip_count['structure']) + '\n'
    analysis_text += '\tHetero:  ' + str(ip_count['hetero']) + '\tHomo:  ' + str(ip_count['homo']) + '\n'
    analysis_text += '\nStructure discovery success (%):  ' + str(sdr) + '\n'
    analysis_text += '\tAmong all trials:\n'
    analysis_text += '\t\tHetero:  ' + str(sdr_hetero_t) + '\tHomo:  ' + str(sdr_homo_t) + '\n'
    analysis_text += '\tAmong individual trials:\n'
    analysis_text += '\t\tHetero:  ' + str(sdr_hetero_i) + '\tHomo:  ' + str(sdr_homo_i) + '\n'
    analysis_text += '\nTotal number of structures discovered:  ' + str(total_num_of_structures) + '\n'
    analysis_text += '\tHetero:  ' + str(hetero_structures) + '\tHomo:  ' + str(homo_structures) + '\n'
    analysis_text += '\nStructures discovered per trial (avg):  ' + str(sdpt) + '\n'
    analysis_text += '\tAmong all trials of same type:\n'
    analysis_text += '\t\tHetero:  ' + str(sdpt_hetero_t) + '\tHomo:  ' + str(sdpt_homo_t) + '\n'
    analysis_text += '\tAmong successful trials of same type:\n'
    analysis_text += '\t\tHetero:  ' + str(sdpt_hetero_i) + '\tHomo:  ' + str(sdpt_homo_i) + '\n\n'

    analysis_path = os.path.join(results_dir, 'ipmof_summary.txt')
    with open(analysis_path, 'w') as a:
        a.write(analysis_text)

        if table:
            header = ['MOF Combination', 'Num of Structures', 'Min Energy', 'Collision', 'Rotation']
            if sortby == 'structure':
                sort_key = 1
            elif sortby == 'energy':
                sort_key = 2
            elif sortby == 'collision':
                sort_key = 3
            sorted_table = sorted(table_lines, key=lambda x: x[1], reverse=True)
            a.write(tabulate(sorted_table, headers=header))
