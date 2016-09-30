#!/bin/bash

#PBS -j oe
#PBS -N ipmof
#PBS -q ishared
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -l mem=1GB
#PBS -S /bin/bash

# accepts a parameter stay_alive if you don't want the worker to exit immediately after all jobs
# are complete
# use like `qsub -v stay_alive=1`

echo JOB_ID: $PBS_JOBID JOB_NAME: $PBS_JOBNAME HOSTNAME: $PBS_O_HOST
echo start_time: `date`

## Load python virtual environment
. /user_path/venv/ipmof/bin/activate

cd $PBS_O_WORKDIR
if [ -z "$stay_alive" ]; then
  options=''
else
  options="--stay-alive"
fi
sjs_launch_workers $PBS_NUM_PPN $options

# workaround for .out / .err files not always being copied back to $PBS_O_WORKDIR
cp /var/spool/torque/spool/$PBS_JOBID.OU $PBS_O_WORKDIR/$PBS_JOBID$(hostname)_$$.out
cp /var/spool/torque/spool/$PBS_JOBID.ER $PBS_O_WORKDIR/$PBS_JOBID$(hostname)_$$.err

exit
