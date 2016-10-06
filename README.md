# IPMOF
### Discovering interpenetration in MOFs

## Installation
IPMOF uses Python 3.5.1 with required libraries listed in requirements.txt file.

To install dependencies just type:

`pip install -r requirements.txt`

## Usage
IPMOF reads structure files from the mof folder in root directory.

The default file reader is ASE and IPMOF supports all the file formats supported by ASE.
A list of the formats supported by ASE is given here: https://wiki.fysik.dtu.dk/ase/ase/io/io.html

In order to test IPMOF there is one MOF file (MOF-5 / REFCODE: SAHYIK) in the 'mof' directory.

## Energy Map
IPMOF is separated into two parts where first a 3D potential energy map is constructed and then
interpenetration is tested making use of this energy map. This allows user to rapidly test
interpenetration many times with different simulation parameters or different MOF combinations
using the same energy map.

Energy map is exported as a numpy array consisting of 3D coordinates and energy values.
The energy values correspond to energy penalties for inserting an atom to a particular position.
As a result energy map is dependent on types of atoms. There are four types of atom lists that
can be used to generate the atom list for a particular MOF:
- 'full': Full atom list in the force field parameters database. (103 atoms)
- 'uniq': Unique atoms for a given list of MOFs
- 'dummy': Single dummy atom with predefined force field parameters
- 'qnd': Simplified atom list consisting of 1 dummy atom and 10 most frequent atoms in MOFs

The runtime of the energy map and interpenetration codes depend on this parameter in the following
fashion: 'full' > 'uniq' > 'qnd' > 'dummy'. The accuracy of the results depend on the set of MOFs
however in general it follows the same trend: 'full' ~ 'uniq' > 'qnd' > 'dummy'. All list types
except for 'uniq' allows checking interpenetration of the energy map MOF with any other MOF.
The 'uniq' energy map uses only the atoms in the given MOF set, therefore any other MOF that
has other types of atoms requires new energy map. 'uniq' and 'qnd' enery map types were designed
to work with large scale screening calculations to have one type of energy map for large
amount of MOF combinations. 'full' energy map can also be used in the same application to
increase accuracy however that would result in slower runtime and bigger file sizes for energy maps.

To generate energy map type following in a terminal window:

`python ipmof_energymap.py`

By default this will create energy maps for each MOF file in ~/mof directory.
The atom list and energy map are stored as a numpy array. This can be changed to a human
readable format (yaml) by changing the 'energy_map_type' simulation parameter to 'yaml'.

## Interpenetration
After energy map is generated interpenetration of two MOFs can be tested.
Energy map is required only for one of the MOFs (should be in ~/energymap) and structure files
(such as 'cif') are required for both MOFs (should be in ~/mof).

To test intepenetration type following in a terminal window:

`python ipmof_interpenetration.py`

Simulation summary, simulation parameters, and information on discovered structures will be
exported to ~/results/*Structure1_Structure2*/results.yaml. Methods to analyze results are
included in ~ipmof/analysis.py library. The structures discovered will be exported to
~/results/*Structure1_Structure2*.

### Simulation Parameters
Default simulation parameters are read from ~/ipmof/parameters.py.

To change the parameters either modify the parameters in that file or create a 'yaml'
file using the functions in ipmof.parameters library. Necessary functions for creating
simulation parameters input files and reading them are given in this library with instructions.
Also, a sample simulation parameters file (sim_par_sample.yaml) is given in ~/settings.
This file can also be used by modifying its name to sim_par.yaml in the same directory.

### Simulation directories
Default directories to read input files such as energy maps, MOF files and output result files
are read from ~/ipmof/parameters.py.

If you wish to use different directories export sim_dir.yaml file into ~/setting folder with
the directories you want to use. You can use methods in ~/ipmof/parameters.py to export this file.

### Documents
To get a better understanding of the functions used in ipmof libraries and/or analyze your
results you can go through jupyter notebook files in '/doc' directory.

This directory also contains force field parameters (UFF and DRE) and supplementary information for
MOFs in CoRE database (10.1021/cm502594j).

### Contact
For any questions, ideas or feedback please contact me!

mail: kbs37@pitt.edu

http://wilmerlab.com/
