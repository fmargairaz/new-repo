
USEFUL COMMAND ON KINGSPEAK FOR INTEL COMPILER
----------------------------------------------
+ old system
source /uufs/chpc.utah.edu/sys/pkg/intel/ics/bin/compilervars.sh intel64
source /uufs/kingspeak.peaks/sys/pkg/mvapich2/std_intel/etc/mvapich2.sh

+ new using module...
module load intel


EXECUTION COMMAND ON KINGSPEAK
------------------------------
+ interactive job
srun --time=00:10:00 --nodes=1 --account=owner-guest --partition=kingspeak-guest --pty /bin/bash -l

+ job and queue checks
squeue -u $USER
scontrol show job

+ kingspeak partition
 - general alocation queue 
#SBATCH --account=calaf
#SBATCH --partition=kingspeak

 - guest node for debug and tests
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest

 - private nodes
#SBATCH --account=calaf-kp
#SBATCH --partition=calaf-kp