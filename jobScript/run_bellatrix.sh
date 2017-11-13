#!/bin/bash
#
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1 
#PBS -l walltime=05:00:00
#PBS -j oe
#PBS -W group_list=seses
#PBS -N fft_bench_test

LD_LIBRARY_PATH=/home/giometto/fft_pjt_fmargairaz/hdf5/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

I_MPI_PIN_DOMAIN=socket
export I_MPI_PIN_DOMAIN

ulimit -s unlimited
ulimit -c unlimited

cd /home/giometto/fft_pjt_fmargairaz/code_fft_test_mpi

echo ' ****** START OF JOB ******'
mpirun ./test_code_mpi.x < inParam_node1_mpi1
echo ' ****** END OF JOB ****** '


