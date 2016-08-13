# Create simulation input files for IPMOF simulations
import os

import yaml

# Simulation Parameters Data:
sim_par_data = {'structure_energy_limit': 3E8,   # Maximum allowed potential energy for structure
                'atom_energy_limit': 3E6,        # Maximum allowed potential energy for atom
                'rotation_limit': 20,            # Total number of rotations for each point
                'rotation_freedom': 30,          # Increments of rotation (degrees)
                'summary_percent': 5,            # Percentage increment to acquire summary data
                'cut_off': 12,                   # Cut-off radius for interpenetration (Angstrom)
                'ext_cut_off': 50,               # Cut-off radius for checking extension (Angstrom)
                'grid_size': 1,                  # Grid size for potential energy map (Angstrom)
                'force_field': 'uff',            # Force field selection for LJ ('uff' or 'dre')
                'core_database': False,          # Use CoRE database information or not
                'energy_map_atom_list': 'qnd',   # Atom list for energy map ('full', 'uniq', 'dummy', 'qnd')
                'energy_map_type': 'numpy',      # Energy map file format ('numpy' or 'yaml')
                'export_structures': 1,          # Number of min. energy structures to export
                'export_format': 'xyz',          # Export structure file format
                'export_pbc': True,              # Export coordinates after applying PBC
                'export_single': True,           # Export structures with original atom names
                'export_single_color': True,     # Export structures with each I.P. layer colored
                'export_packed': True,           # Export structures with packed unit cell
                'export_packed_color': True,     # Export structures with packed unit cell and color
                'export_summary': True           # Export summary information for interpenetration
                }

# Working Directories:
main_dir = os.getcwd()
python_lib_dir = os.path.join(main_dir, 'ipmof')
force_field_path = os.path.join(main_dir, 'doc', 'FF_Parameters.xlsx')
core_path = os.path.join(main_dir, 'doc', 'CoRE.xlsx')
core_mof_dir = r'CoRE MOFs database direcory'
mof_dir = os.path.join(main_dir, 'mof')
energy_map_dir = os.path.join(main_dir, 'energymap')
export_dir = os.path.join(main_dir, 'results')
if not os.path.isdir(export_dir):
    os.mkdir(export_dir)

# Simulation Directories Data:
sim_dir_data = {'main_dir': main_dir,
                'python_lib_dir': python_lib_dir,
                'force_field_path': force_field_path,
                'core_path': core_path,
                'core_mof_dir': core_mof_dir,
                'mof_dir': mof_dir,
                'energy_map_dir': energy_map_dir,
                'export_dir': export_dir,
                }


def export_sim_par(inp_dir=main_dir):
    """
    Export simulation parameters input file to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_par_path = os.path.join(input_dir, 'sim_par.yaml')
    with open(sim_par_path, 'w') as sim_par_file:
        yaml.dump(sim_par_data, sim_par_file, default_flow_style=False)


def export_sim_dir(inp_dir=main_dir):
    """
    Export simulation directories to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_dir_path = os.path.join(input_dir, 'sim_dir_linux.yaml')
    with open(sim_dir_path, 'w') as sim_dir_file:
        yaml.dump(sim_dir_data, sim_dir_file, default_flow_style=False)

# Get directories for simulation parameters and directories files
# sim_par_path = os.path.join(main_dir, 'sim_par.yaml')
# sim_dir_path = os.path.join(main_dir, 'sim_dir_linux.yaml')

# Read sim par yaml file
# sim_par = yaml.load(open(sim_par_path, 'r'))
# sim_dir = yaml.load(open(sim_dir_path, 'r'))


def export_init_txt(mof_list, sim_par=sim_par_data, sim_dir=sim_dir_data):
    """
    Export init.txt file to results folder.
    Contains information about simulation parameters and selected structures.
    """
    init_text = '--------------------- IPMOF ---------------------\n'
    init_text += '------------- SIMULATION PARAMETERS -------------\n'
    max_len = len('structure_energy_limit') + 3
    for par in sim_par:
        len_par_name = len(str(par))
        empty_space = ' ' * (max_len - len_par_name)
        init_text += par + ':' + empty_space + str(sim_par[par]) + '\n'
    init_text += '-------------------------------------------------\n'

    init_text += 'Starting interpenetration with a total of ' + str(len(mof_list)) + ' MOFs:\n'
    for m_i, m in enumerate(mof_list):
        init_text += str(m_i + 1) + '\t' + str(m) + '\n'
    init_text += '-------------------------------------------------\n'

    init_file_path = os.path.join(sim_dir['export_dir'], 'init.txt')
    if os.path.exists(init_file_path):
        os.remove(init_file_path)

    with open(init_file_path, 'w') as init_file:
        init_file.write(init_text)


def export_summary_txt(export_dir, summary, base_mof, mobile_mof):
    """
    Export summary.txt file to specific interpenetration folder in results folder.
    Contains information about interpenetration simulation run.
    """
    summary_text = '--------------- IPMOF SUMMARY ---------------\n'
    summary_text += 'Base -> ' + base_mof.name + '\nMobile -> ' + mobile_mof.name + '\n'
    summary_text += '---------------------------------------------\n'
    summary_text += 'Percent \t Structure \t Trial \n'
    for summary_index, percent in enumerate(summary['percent']):
        summary_text += str(percent) + ' \t\t ' + str(summary['structure_count'][summary_index])
        summary_text += ' \t\t ' + str(summary['trial_count'][summary_index]) + ' \n'
    summary_text += '---------------------------------------------\n'

    summary_file_path = os.path.join(export_dir, 'summary.txt')
    if os.path.exists(summary_file_path):
        os.remove(summary_file_path)
    with open(summary_file_path, 'w') as summary_file:
        summary_file.write(summary_text)
