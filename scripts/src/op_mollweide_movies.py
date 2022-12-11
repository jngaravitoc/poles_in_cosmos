"""
Master script to analyze the FIRE simulations for the orbital poles project
(github/jngaravitoc/poles_in_cosmos)


This script has been tested with sims: m12b, m12i

Main functionalities:
   - Make plots
    - Density plots of the DM and stellar distribution in several projections 
    - Mollweide plots of particles and subhalos positions in Galactocentric
      coordinates.
    - Mollweide plots of the orbital poles
   - Perform analysis
    - Correlation function analysis

Dependencies:
  - scipy
  - numpy 
  - Gizmo Analysis
  - Halo tools
  - pynbody
  - Astropy
  - nba 

Author: Nico Garavito-Camargo
Github: jngaravitoc

TODO:
- Remove satellite subhalos

"""

#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import sys
import pynbody

sys.path.append("/mnt/home/ecunningham/python")
plt.style.use('~/matplotlib.mplstyle')
import gizmo_analysis as ga
import halo_analysis as halo
import nba

# 
import pynbody_routines  as pr 
import plotting as pl
from io_gizmo_pynbody  import FIRE
import analysis as an

if __name__ == "__main__":
    
    
    snap_init = 450
    snap_final = 451
    sim='m12f'
    sats = False
    rmin = 50
    rmax = 300
    
    snap_times = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt'.format(sim)
    times = np.loadtxt(snap_times, usecols=3)
    
    
    m12 = FIRE(sim, remove_satellite=True)
    for k in range(snap_init, snap_final):
        # face on particle data halo
        hfaceon = m12.rotated_halo(k)

        # Satellite orbit
        subhalos_faceon, satellite_faceon = m12.subhalos_rotated(k)

        
        pos_dm = hfaceon.dark['pos']
        f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
        vel_dm = hfaceon.dark['vel']*f
        dist_dm = np.sqrt(np.sum(pos_dm**2, axis=1))

        dist_cut1 = np.where((dist_dm> rmin) & (dist_dm< rmax)) 
        #ist_cut2 = np.where((dist_dm> 300) & (dist_dm< 600)) 
        
        dm_kinematics1 = nba.kinematics.Kinematics(pos_dm[dist_cut1],  vel_dm[dist_cut1])
        #dm_kinematics2 = nba.kinematics.Kinematics(pos_dm[dist_cut2],  vel_dm[dist_cut2])
        
        vel_galactic = dm_kinematics1.vel_cartesian_to_galactic()
        
        vis = nba.Visuals(pos_dm[dist_cut1])
        vis.all_sky_galactic(100, 50, vel_galactic[0])

        dm_figname1 = "../plots/exploration/{}_OP_satellite_faceon_{:03d}.png".format(sim, k)
        fig_title = "{} satellite poles DM; {}-{} kpc; t={:.2f}  Gyr".format(sim, 50, 300, t_snap[k] )

        pl.mollweide_projection(dm_OP_l_host2, dm_OP_b_host2, sat_OP_l_host[k-300], sat_OP_b_host[k-300], title=title_dm2, bmin=0, bmax=500, nside=40, smooth=1, figname=dm_figname2)

        #sat_kinematics = nba.kinematics.Kinematics(satellite_faceon.dark['pos'],  satellite_faceon.dark['vel'])
        #sat_l_host, sat_b_host = sat_kinematics.pos_cartesian_to_galactic()
        #sat_OP_l_host, sat_OP_b_host = sat_kinematics.orbpole()

        #dm_figname1 = "../plots/exploration/{}_galacitc_faceon_{:03d}.png".format(sim, k)
                             
        #title_dm1 = "{} satellite poles DM; {}-{} kpc; t={:.2f}  Gyr".format(sim, 50, 300, t_snap[k] )
        #title_dm2 = "{} satellite poles DM; {}-{} kpc; t={:.2f}  Gyr".format(sim, 300, 600, t_snap[k] )
        #dm_figname1 =  "../plots/exploration/{}_OP_satellite_faceon_{:03d}_50_300.png".format(sim, k)
        #dm_figname2 =  "../plots/exploration/{}_OP_satellite_faceon_{:03d}_300_600.png".format(sim, k)

        #pl.mollweide_projection(dm_OP_l_host1, dm_OP_b_host1, sat_OP_l_host[k-300], sat_OP_b_host[k-300], 
         #                    title=title_dm1, bmin=100, bmax=1000,
         #                    nside=40, smooth=1, figname=dm_figname1)
     
        #pl.mollweide_projection(dm_OP_l_host2, dm_OP_b_host2, sat_OP_l_host[k-300], sat_OP_b_host[k-300], 
        #                     title=title_dm2, bmin=0, bmax=500,
        #                     nside=40, smooth=1, figname=dm_figname2)
