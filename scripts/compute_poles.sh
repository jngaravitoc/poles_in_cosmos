#!/bin/bah

module load python 
ncores=24
sim=m12b
figname=OP_dark_particles_shells_no_rot_no_substructure_
rmin=50
rmax=500
snapi=500
snapf=601

python FIRE_analysis.py --ncores=$ncores  --sim=$sim --figname=$figname --rmin=$rmin --rmax=$rmax --i=$snapi --f=$snapf
