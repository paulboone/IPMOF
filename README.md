# IPMOF
## Discovering interpenetration in MOFs

## Installation
IPMOF uses Python 3.5.1 with required libraries listed in requirements.txt file.

To install dependencies just type:

`pip install -r requirements.txt`

If you do not wish to install all the libraries in the requirements file 'ase' and 'openbabel'
libraries can be skipped. However only '.mol2' format is supported without these libraries.
It is recommended to install 'ase' library for multiple format recognition.

## Usage
IPMOF reads structure files from the mof folder in root directory.

The default file reader is ASE and IPMOF supports all the file formats supported by ASE.
A list of the formats supported by ASE is given here: https://wiki.fysik.dtu.dk/ase/ase/io/io.html

IPMOF also has a built-in '.mol2' reader which does not requiee additional libraries.

In order to test IPMOF there is one MOF file (MOF-5 / REFCODE: SAHYIK) in the 'mof' directory.

### Energy Map
IPMOF is separated into two parts where first a 3D potential energy map is constructed and
interpenetration is tested making use of this energy map. This allows user to rapidly test
interpenetration many times with different simulation parameters or different MOF combinations
by using the same energy map.

Energy map is exported as a numpy array consisting of 3D coordinates and energy values.
The energy values correspond to energy penalties for inserting an atom to a particular position.
As a result energy map is dependent on types of atoms. There are four types of atom lists that
can be used to generate the atom list for a particular MOF. All list types except for 'uniq' allows
checking interpenetration of the energy map MOF with any other MOF.

To generate energy map type following in a terminal window:

`python ipmof_energymap.py`

By default this will create energy maps for each MOF file in ~/mof directory.
The atom list and energy map are stored as a numpy array. This can be changed to a human
readable format (yaml) by changing the 'energy_map_type' simulation parameter to 'yaml'.

### Interpenetration
After energy map is generated interpenetration of two MOFs, where is always the MOF
the energy map is created for, can be tested.

To test intepenetration type following in a terminal window:

`python ipmof_interpenetration.py`

The simulation parameters used in the simulation are exported to ~/results/init.txt
Simulation summary and structures discovered will be exported to ~/results/*Structure1_Structure2*

If you do not wish to export summary file you can change summary_export to False in parameters.
Also you can choose structure export type and format in parameters as well.

### Simulation Parameters
Default simulation parameters are read from ~/ipmof/parameters.py.

To change the parameters either modify the parameters in that file or create a 'yaml'

file using the functions in ipmof.parameters library. Necessary functions for creating

simulation parameters input files and reading them are give in this library with instructions.

### IPMOF Libraries
To get a better understanding of the functions used in ipmof libraries you can go through
jupyter notebook files in '/doc' directory.

A notebook version of IPMOF is also present in this directory which can be used to make
changes in the main IPMOF algorithm.
