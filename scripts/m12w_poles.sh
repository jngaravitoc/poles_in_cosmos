#!/bin/bash 

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ngaravito@flatironinstitute.org
#SBATCH --time=3:00:00
#SBATCH --job-name m12w
#SBATCH -N1 --ntasks-per-node=20 
#SBATCH -e stderr_m12w.txt
#SBATCH -o stdout_m12w.txt


module purge
module load python

echo
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo


ncores=20
sim=m12w
snapi=300
snapf=600
partype=1 # 1) DM 2) stars 3) both
#figname=OP_disk_no_rot_rleq_50


figname=OP_dark_no_rot_r50_150
rmin=50
rmax=150
outpath="/mnt/home/nico/projects/poles_in_cosmos/plots/${sim}/OP_50_150_no_subs/"

python FIRE_analysis.py --ncores=$ncores  --sim=$sim --figname=$figname --outpath=$outpath  --rmin=$rmin --rmax=$rmax --i=$snapi --f=$snapf --remove_subs=True --partype=$partype

# Make movie
sh movies.sh -i "${outpath}${sim}_${figname}_" -o "${outpath}${sim}${figname}.mp4" -s $snapi

# ----------------------

figname=OP_dark_no_rot_r150_300
rmin=150
rmax=300
outpath="/mnt/home/nico/projects/poles_in_cosmos/plots/${sim}/OP_150_300_no_subs/"

python FIRE_analysis.py --ncores=$ncores  --sim=$sim --figname=$figname --outpath=$outpath  --rmin=$rmin --rmax=$rmax --i=$snapi --f=$snapf --remove_subs=True --partype=$partype

# Make movie
sh movies.sh -i "${outpath}${sim}_${figname}_" -o "${outpath}${sim}${figname}.mp4" -s $snapi

# ----------------------

figname=OP_dark_no_rot_r300_500
rmin=300
rmax=500
outpath="/mnt/home/nico/projects/poles_in_cosmos/plots/${sim}/OP_300_500_no_subs/"

python FIRE_analysis.py --ncores=$ncores  --sim=$sim --figname=$figname --outpath=$outpath  --rmin=$rmin --rmax=$rmax --i=$snapi --f=$snapf --remove_subs=True --partype=$partype

# Make movie
sh movies.sh -i "${outpath}${sim}_${figname}_" -o "${outpath}${sim}${figname}.mp4" -s $snapi


# ----------------------

figname=OP_dark_no_rot_r50_500
rmin=50
rmax=300
outpath="/mnt/home/nico/projects/poles_in_cosmos/plots/${sim}/OP_50_500/"

python FIRE_analysis.py --ncores=$ncores  --sim=$sim --figname=$figname --outpath=$outpath  --rmin=$rmin --rmax=$rmax --i=$snapi --f=$snapf --remove_subs=False --partype=$partype

# Make movie
sh movies.sh -i "${outpath}${sim}_${figname}_" -o "${outpath}${sim}${figname}.mp4" -s $snapi



# ----------------------

figname=OP_disk_no_rot_r0_50
rmin=0
rmax=50
outpath="/mnt/home/nico/projects/poles_in_cosmos/plots/${sim}/disk/"
partype=2 # 1) DM 2) stars 3) both

python FIRE_analysis.py --ncores=$ncores  --sim=$sim --figname=$figname --outpath=$outpath  --rmin=$rmin --rmax=$rmax --i=$snapi --f=$snapf --remove_subs=False --partype=$partype

# Make movie
sh movies.sh -i "${outpath}${sim}_${figname}_" -o "${outpath}${sim}${figname}.mp4" -s $snapi
