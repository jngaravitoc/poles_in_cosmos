#!/bin/bash 

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ngaravito@flatironinstitute.org
#SBATCH --time=8:00:00
#SBATCH --job-name m12b_poles
#SBATCH -N1 --ntasks-per-node=1
#SBATCH -e stderr.txt
#SBATCH -o stdout.txt

module purge
module load slurm
module load python


echo
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo

#python  /mnt/home/nico/projects/poles_in_cosmos/scripts/m12b_orbital_poles_time_evolution.py 
python /mnt/home/nico/projects/poles_in_cosmos/scripts/gizmo_pynbody_analysis_workflow.py
#python  /mnt/home/nico/projects/poles_in_cosmos/scripts/mwlmc5_orbital_poles_computation.py 

##mpirun -np  $SLURM_NPROCS  /mnt/home/nico/codes/gadget4/Gadget4_grav_MO3 /mnt/home/nico/gadget_runs/isolated_halo/param_ics2_grav_MO3.txt 2 223 
