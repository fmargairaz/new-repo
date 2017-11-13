#!/bin/bash

#SBATCH --time=65:00:00 
#SBATCH --nodes=4 
#SBATCH -o slurm-%j.out-%N 
#SBATCH --ntasks=16
#SBATCH --account=calaf
#SBATCH --partition=kingspeak
#SBATCH -J JOBNAME
#SBATCH --mail-user=USER-EMAIL
#SBATCH --mail-type=FAIL,BEGIN,END 

# Prologue
echo '************************ PROLOGUE ************************'
echo 'Hi, jobID: '$SLURM_JOBID
date
echo '----------------------------------------------------------'
echo 'setting environment'

# load module

I_MPI_PIN_DOMAIN=socket
export I_MPI_PIN_DOMAIN
OMP_NUM_THREADS=4
export OMP_NUM_THREADS

echo '----------------------------------------------------------'

cd /scratch/kingspeak/serial/....

echo ' ****** START OF JOB ******'
mpirun -np $SLURM_NTASKS ./les_code_mpi.x > std_output
echo ' ****** END OF JOB ****** '
