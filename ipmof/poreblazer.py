# PoreBlazer functions for importing and exporting PoreBlazer simulation files.
# Sarkisov et al. Molecular Simulation, 37, 1248-1257 (2011) - DOI:10.1080/08927022.2011.592832
# Date: June 2016
# Author: Kutay B. Sezginel
import os
import shutil
import matplotlib.pyplot as plt


class PB:
    """
    PB class for importing and exporting PoreBlazer (DOI:10.1080/08927022.2011.592832)
    simulation files.
    """
    def source(self, source_dir):
        """
        Initialize source directory for PoreBlazer simulation files.
        """
        self.source_dir = source_dir

    def input(self, structure_name, uc_size, uc_angle, export_dir):
        """
        Creates input.dat file according to unit cell information of given structure.
        """
        pb_input_dir = os.path.join(export_dir, 'input.dat')
        with open(pb_input_dir, 'w') as pb_input_file:
            pb_input_file.write(structure_name + '\n')

            uc_size_line = str(uc_size[0]) + ' ' + str(uc_size[1]) + ' ' + str(uc_size[2]) + '\n'
            pb_input_file.write(uc_size_line)

            uc_angle_line = str(uc_angle[0]) + ' ' + str(uc_angle[1]) + ' ' + str(uc_angle[2]) + '\n'
            pb_input_file.write(uc_angle_line)

    def initialize(self, export_dir):
        """
        Copies all the files in the source directory to export directory.
        """
        source_dir = self.source_dir
        for inp_file in os.listdir(source_dir):
            inp_dir = os.path.join(source_dir, inp_file)
            shutil.copy(inp_dir, export_dir)

    def JobSH(self, job_dir, job_name, queue, wall_time):
        """
        Creates job submission file to run PoreBlazer simulations on Frank.
        """
        jobsh = []
        jobsh.append('#!/bin/bash')
        jobsh.append(' ')
        jobsh.append('#PBS -j oe')
        jobsh.append('#PBS -N ' + job_name)
        jobsh.append('#PBS -q ' + queue)
        jobsh.append('#PBS -l nodes=1:ppn=1')
        jobsh.append('#PBS -l walltime=' + wall_time)
        jobsh.append('#PBS -S /bin/bash')
        jobsh.append(' ')
        jobsh.append('echo JOB_ID: $PBS_JOBID JOB_NAME: $PBS_JOBNAME HOSTNAME: $PBS_O_HOST')
        jobsh.append('echo "start_time: `date`"')
        jobsh.append(' ')
        jobsh.append('cd $PBS_O_WORKDIR')
        jobsh.append('./ihome/cwilmer/kbs37/PoreBlazer/poreblazer.exe < input.dat > results.txt')
        jobsh.append('echo end_time: `date`')
        jobsh.append(' ')
        jobsh.append('cp /var/spool/torque/spool/$PBS_JOBID.OU $PBS_O_WORKDIR/$PBS_JOBID.out')
        jobsh.append('cp /var/spool/torque/spool/$PBS_JOBID.ER $PBS_O_WORKDIR/$PBS_JOBID.err')
        jobsh.append(' ')
        jobsh.append('exit')
        jobsh_dir = os.path.join(job_dir, 'job.sh')
        with open(jobsh_dir, 'w') as jobsh_file:
            for line in jobsh:
                jobsh_file.write(line + '\n')

    def read_results(self, pb_results_dir):
        """
        Reads results.txt file from PoreBlazer simulation output.
        Collects surface area (m2/g), density(g/cm3), and geometric pore volume (cm3/g).

        Example usage:
         >>> pb_res = PB.read_results(pb_res_dir)
         ... pb_res = {'sa': *surface_area, 'pv': *pore_volume, 'ro': *crystal_density}
        """
        pb_results_dir = os.path.join(pb_results_dir, 'results.txt')
        with open(pb_results_dir, 'r') as pb_results_file:
            pb_results_lines = pb_results_file.readlines()
        for line in pb_results_lines:
            if 'surface area per mass' in line:
                sa = float(line.split()[-1])
            if 'System density' in line:
                ro = float(line.split()[-1])
            if '(point accessible) volume in cm^3/g' in line:
                pv = float(line.split()[-1])
        pb_res = {}
        pb_res['sa'] = sa
        pb_res['ro'] = ro
        pb_res['pv'] = pv

        return pb_res

    def read_psd(self, psd_results_dir):
        """
        Reads psd.txt file from PoreBlazer simulation output.
        Collects pore size and frequency data from file and calculates dominant pore diameter.
        Output dictionary can be plotted with plot_psd function.

        Example usage:
         >>> psd_res = PB.read_psd(psd_res_dir)
         >>> psd_res = {'pore_size': [pore size],
                        'frequency': [frequency],
                        'dpd': *dominant pore diameter}
        """
        psd_dir = os.path.join(psd_results_dir, 'psd.txt')
        with open(psd_dir, 'r') as psd_file:
            psd_lines = psd_file.readlines()

        pore_size = []
        frequency = []
        for i, line in enumerate(psd_lines):
            if i > 0:
                pore_size.append(float(line.split()[0]))
                frequency.append(float(line.split()[1]))

        psd = {}
        psd['pore_size'] = pore_size
        psd['frequency'] = frequency
        dpd_index = psd['frequency'].index(max(psd['frequency']))
        psd['dpd'] = psd['pore_size'][dpd_index]

        return psd

    def plot_psd(self, psd):
        """
        Plots pore size distribution results read from output of PoreBlazer simulations.
        """
        import matplotlib.pyplot as plt
        plt.subplots(figsize=(8, 5))
        plt.plot(psd['pore_size'], psd['frequency'])
        plt.title('Pore Size Distribution', y=1.04, fontsize=18)
        plt.ylabel('Frequency', fontsize=18)
        plt.xlabel('Pore Size (Angstrom)', fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        text_x_pos = psd['dpd'] * 1.1
        text_y_pos = max(psd['frequency']) * 0.97
        plt.text(text_x_pos, text_y_pos, 'dpd: ' + str(psd['dpd']), fontsize=14)
        plt.show()
