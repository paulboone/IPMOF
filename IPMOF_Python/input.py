# Create simulation input files for IPMOF simulations
import os
import yaml

sim_par_data = {'structure_energy_limit': 3E3,
                'atom_energy_limit': 3E1,
                'rotation_limit': 20,
                'rotation_freedom': 30,
                'summary_percent': 5,
                'cut_off': 12,
                'ext_cut_off': 50,
                'grid_size': 1,
                'force_field': 'uff'
                }

# Enter directories here
sim_dir_data = {'python_lib_dir': r'/home/kutay/Documents/git/IPMOF/IPMOF_Python',
                'excel_file_path': r'/home/kutay/Documents/Research/IPMOF/FF_Parameters.xlsx',
                'mol2_dir': r'/home/kutay/Documents/Research/MOFs/IPMOF_Python/mol2',
                'input_dir': r'/home/kutay/Documents/git/IPMOF',
                'export_dir': r'/home/kutay/Documents/Research/MOFs/IPMOF_Python/export',
                }


def export_sim_par(inp_dir=input_dir):
    """
    Export simulation parameters input file to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_par_path = os.path.join(input_dir, 'sim_par.yaml')
    sim_par_file = open(sim_par_path, 'w')
    yaml.dump(sim_par_data, sim_par_file, default_flow_style=False)
    sim_par_file.close()


def export_sim_dir(inp_dir=input_dir):
    """
    Export simulation directories to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_dir_path = os.path.join(input_dir, 'sim_dir_linux.yaml')
    sim_dir_file = open(sim_dir_path, 'w')
    yaml.dump(sim_dir_data, sim_dir_file, default_flow_style=False)
    sim_dir_file.close()
