#!/bin/bash 

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ngaravito@flatironinstitute.org
#SBATCH --time=48:00:00
#SBATCH --job-name  m12csubhalos_catalogs
#SBATCH -N1 --ntasks-per-node=8 
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


##python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12f 1 
##python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12f 0 

##python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12r 1 
##python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12r 0 

##python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12m 0 

##python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12w 0 
##python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12w 1 

python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12c 1 
#python  /mnt/home/nico/projects/poles_in_cosmos/scripts/make_subhalo_catalogue.py m12b 1

