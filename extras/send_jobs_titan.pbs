#!/bin/bash
#    Begin PBS directives
#PBS -A ast125
#PBS -N cosmo_1024
#PBS -j oe
#PBS -l walltime=0:05:00
#PBS -l nodes=512
#PBS -q debug
#PBS -l gres=atlas1%atlas2        Optional      The job will require both the atlas1 and atlas2 Lustre® filesystems to be online.
#    End PBS directives and begin shell commands
# source $MODULESHOME/init/bash
module switch PrgEnv-pgi/5.2.82 PrgEnv-gnu
module load cray-hdf5
module load fftw
module load cudatoolkit

cd $MEMBERWORK/ast125/cholla
# cat $PBS_NODEFILE | sort | uniq > hosts.$PBS_JOBID
export OMP_NUM_THREADS=16
date
aprun -n 512 -N 1 -d 16 ./cholla tests/3D/Cosmological_hydro_1024.txt
