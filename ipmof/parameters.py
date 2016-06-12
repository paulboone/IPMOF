# Create simulation input files for IPMOF simulations
import os

import yaml

# Simulation Parameters Data:
sim_par_data = {'structure_energy_limit': 3E10,
                'atom_energy_limit': 3E8,
                'rotation_limit': 30,
                'rotation_freedom': 30,
                'summary_percent': 5,
                'cut_off': 12,
                'ext_cut_off': 50,
                'grid_size': 1,
                'force_field': 'uff',
                'export_structures': 1,
                'export_pbc': False,
                'export_colorify': True
                'export_format': 'xyz'
                }

# Working Directories:
main_dir = os.getcwd()
python_lib_dir = os.path.join(main_dir, 'ipmof')
force_field_path = os.path.join(main_dir, 'doc', 'FF_Parameters.xlsx')
core_path = os.path.join(main_dir, 'doc', 'CoRE.xlsx')
mol2_dir = r'/home/kutay/Documents/Research/MOFs/IPMOF_Python/mol2'
# mol2_dir = r'C:\Kutay\MOFs\IPMOF_Python'
export_dir = os.path.join(main_dir, 'results')
if not os.path.isdir(export_dir):
    os.mkdir(export_dir)

# Simulation Directories Data:
sim_dir_data = {'main_dir': main_dir,
                'python_lib_dir': python_lib_dir,
                'force_field_path': force_field_path,
                'core_path': core_path,
                'mol2_dir': mol2_dir,
                'export_dir': export_dir,
                }


def export_sim_par(inp_dir=main_dir):
    """
    Export simulation parameters input file to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_par_path = os.path.join(input_dir, 'sim_par.yaml')
    sim_par_file = open(sim_par_path, 'w')
    yaml.dump(sim_par_data, sim_par_file, default_flow_style=False)
    sim_par_file.close()


def export_sim_dir(inp_dir=main_dir):
    """
    Export simulation directories to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_dir_path = os.path.join(input_dir, 'sim_dir_linux.yaml')
    sim_dir_file = open(sim_dir_path, 'w')
    yaml.dump(sim_dir_data, sim_dir_file, default_flow_style=False)
    sim_dir_file.close()

# Get directories for simulation parameters and directories files
# sim_par_path = os.path.join(main_dir, 'sim_par.yaml')
# sim_dir_path = os.path.join(main_dir, 'sim_dir_linux.yaml')

# Read sim par yaml file
# sim_par = yaml.load(open(sim_par_path, 'r'))
# sim_dir = yaml.load(open(sim_dir_path, 'r'))
