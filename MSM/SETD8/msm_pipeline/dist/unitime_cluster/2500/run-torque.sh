#!/bin/bash
#  Batch script for mpirun job on cbio cluster.
#
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=120:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q batch
#
# nodes: number of nodes
#   ppn: how many cores per node to use
#PBS -l nodes=1:ppn=10
#
#PBS -l mem=50G
#
# export all my environment variables to the job
##PBS -V
#
# job name (default = name of script file)
#PBS -N cluster-dih_noNflank_2500

if [ -n "$PBS_O_WORKDIR" ]; then 
    cd $PBS_O_WORKDIR
fi

# rm -rf pyemma.log dtrajs clustercenters.npy

# Run the simulation with verbose output:
date
export MPLBACKEND="agg"

export OMP_NUM_THREADS=20
echo $OMP_NUM_THREADS
date
time python cluster.py
date
