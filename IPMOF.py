# IPMOF Main Algorithm
import math
import os
import yaml

input_dir = r'/home/kutay/Documents/git/IPMOF'

sim_par_path = os.path.join(input_dir, 'sim_par.yaml')
sim_dir_path = os.path.join(input_dir, 'sim_dir_linux.yaml')

# Read sim par yaml file
sim_par = yaml.load(open(sim_par_path, 'r'))
sim_dir = yaml.load(open(sim_dir_path, 'r'))

# Load interpenetration python libraries
os.chdir(sim_dir['python_lib_dir'])
from forcefield import read_ff_parameters
from interpenetration import run_interpenetration, check_extension
from crystal import MOF, Packing
from energymap import energy_map, get_mof_list, get_uniq_atom_list,
                      energy_map_index, energy_map_atom_index

# Read excel file containing force field information
uff = read_ff_parameters(sim_dir['excel_file_path'], 'uff')

# Create list of MOFs
mol2_list = get_mof_list(mol2_dir, '.mol2')
print(mol2_list)

base_mof_index = 4
mobile_mof_index = 4

# ------------------------------------- Loop Over Here ---------------------------------------------
# --------------------------------------------------------------------------------------------------
# Read mol2 files and initialize MOF objects
base_mof = MOF()
base_mof.mol2_path = os.path.join(mol2_dir, mol2_list[base_mof_index])
base_mof.initialize()
base_mof.initialize_ff(uff)
print('Base MOF selected as: ', base_mof.name)

mobile_mof = MOF()
mobile_mof.mol2_path = os.path.join(mol2_dir, mol2_list[mobile_mof_index])
mobile_mof.initialize()
mobile_mof.initialize_ff(uff)
print('Mobile MOF selected as: ', mobile_mof.name)

cut_off = sim_par['cut_off']

packing_factor = Packing.factor(base_mof.uc_size, cut_off)
uc_vectors = Packing.uc_vectors(base_mof.uc_size, base_mof.uc_angle)
trans_vec = Packing.translation_vectors(packing_factor, uc_vectors)
base_mof.packed_coors = Packing.uc_coors(trans_vec, packing_factor, uc_vectors, base_mof.atom_coors)
base_mof.edge_points = Packing.edge_points(uc_vectors)

print('Base MOF unit cell: ', base_mof.uc_size)
print('Packing factor:', base_mof.packing_factor)
print('Num of coor :', len(base_mof.packed_coors)*len(base_mof.packed_coors[0]))

atom_list = get_uniq_atom_list([mobile_mof])
print('Calculating emap for', base_mof.name, 'with atoms:', atom_list['atom'])

# Calculate energy map
emap = energy_map(base_mof, atom_list, cut_off, 1)

# Run interpenetration
summary, new_structures = run_interpenetration(sim_par, base_mof, mobile_mof, emap, atom_list)

# Check for extension
collision = check_extension(sim_par, base_mof, mobile_mof, emap, atom_list, rotation_info)

if not collision:
    # Export structure

# --------------------------------------------------------------------------------------------------
