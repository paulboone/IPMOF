# IPMOF
## Discovering interpenetration in MOFs

### Installation
IPMOF uses Python 3.5.1 with required libraries listed in requirements.txt file.

To install dependencies just type:

`pip install -r requirements.txt`

### Usage
IPMOF reads structure files from the mof folder in root directory.

The default format for structure files is '.mol2'.

In order to test IPMOF there is one MOF file (MOF-5) in the 'mof' directory.

To run simulations in your terminal window just type:
```
python IPMOF.py or
python IPMOF.py > results.txt
```
if you want to export output into a file.

This will simply run interpenetration simulation for MOF-5.
The simulation parameters used in the simulation are exported to './results/init.txt'

Discovered structure coordinates are also exported to
'/results/*Structure1_Structure2*/*Structure1_Structure2.xyz*'.

Simulation summary for each interpenetration is exported to
'/results/*Structure1_Structure2*/summary.txt'

If you do not wish to export summary file you can change summary_export to False in parameters.

The sample simulation (MOF-5 + MOF-5) takes around 18 minutes
on Ubuntu 16.04 LTS with Intel® Core™ i3-6100U CPU @ 2.30GHz × 4.

#### Simulation Parameters
Default simulation parameters are read from '/ipmof/parameters.py'.

To change the parameters either modify the parameters in that file or create a 'yaml'

file using the functions in ipmof.parameters library. Necessary functions for creating

simulation parameters input files and reading them are give in this library with instructions.

#### IPMOF Libraries
To get a better understanding of the functions used in ipmof libraries you can go through
jupyter notebook files in '/doc' directory.

A notebook version of IPMOF is also present in this directory which can be used to make
changes in the main IPMOF algorithm.
