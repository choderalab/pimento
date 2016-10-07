#!/bin/tcsh
#  Batch script for MPI GPU job on the cbio cluster
#  utilizing 4 GPUs, with one thread/GPU
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=72:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify GPU queue
#PBS -q gpu
#
# nodes: number of nodes
#   ppn: number of processes per node
#  gpus: number of gpus per node
#  GPUs are in 'exclusive' mode by default, but 'shared' keyword sets them to shared mode.
#PBS -l nodes=1:ppn=1:gpus=1
#
# export all my environment variables to the job
#PBS -V
#
# job name (default = name of script file)
#PBS -N SET8_SAM
#
# mail settings (one or more characters)
# email is sent to local user, unless another email address is specified with PBS -M option 
# n: do not send mail
# a: send mail if job is aborted
# b: send mail when job begins execution
# e: send mail when job terminates
#PBS -m n
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
##PBS -o myoutput

# Change to working directory used for job submission
cd $PBS_O_WORKDIR

# Launch GPU job.
python -u equilibrate_package_SET8_SAM.py > equilibrate_package_SET8_SAM.log
