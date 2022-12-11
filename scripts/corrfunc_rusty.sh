#!/bin/bash 

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ngaravito@flatironinstitute.org
#SBATCH --time=24:00:00
#SBATCH --job-name 50_300_corrfunc_dm_m12r
#SBATCH --array=0-30 -N1 --ntasks-per-node=1 
#SBATCH -e stderr.%j.%A.%a.%N.txt
#SBATCH -o stdout.%j.%A.%a.%N.txt

module purge
module load slurm
module load python


echo
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo

START_NUM=$(( ($SLURM_ARRAY_TASK_ID * 10) + 300))
END_NUM=$(($START_NUM+10))

python  /mnt/home/nico/projects/poles_in_cosmos/scripts/src/correlation_functions_particles.py $START_NUM $END_NUM $SLURM_CPUS_ON_NODE


