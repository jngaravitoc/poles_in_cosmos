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
#plt.style.use('~/matplotlib.mplstyle')
import gizmo_analysis as ga
import halo_analysis as halo
import nba

# 
import pynbody_routines  as pr 
import plotting as pl
from io_gizmo_pynbody  import FIRE
import analysis as an

if __name__ == "__main__":
    
    
    snap_init = 352
    snap_final = 354
    sim='m12b'
    #sats = False
    rmin = 50
    rmax = 300
    ptype = ['star']

    snap_times = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt'.format(sim)
    times = np.loadtxt(snap_times, usecols=3)
    plot_type = ['cartesian_projection'] # vr_mollweide, orbital_poles 
    
    m12 = FIRE(sim, remove_satellite=True, rm_stellar_sat=True)

    for k in range(snap_init, snap_final,2):
        # face on particle data halo
        hfaceon = m12.rotated_halo(k)

        # Satellite orbit
        subhalos_faceon, satellite_faceon = m12.subhalos_rotated(k)

        f = 1* (u.km/u.s).to(u.kpc/u.Gyr)

        if 'star' in ptype:
            pos_s = hfaceon.star['pos']
            vel_s = hfaceon.star['vel']*f
            dist_s = np.sqrt(np.sum(pos_s**2, axis=1))
            print(len(pos_s))
        if 'dark' in ptype:
            pos_dm = hfaceon.dark['pos']
            vel_dm = hfaceon.dark['vel']*f
            dist_dm = np.sqrt(np.sum(pos_dm**2, axis=1))
        if ("vr_mollweide" or  "orbital_poles" or "rho_mollweide") in plot_type:
            kinematics1 = nba.kinematics.Kinematics(pos_dm[dist_cut1],  vel_dm[dist_cut1])
            #kinematics2 = nba.kinematics.Kinematics(pos[dist_cut2],  vel[dist_cut2])
               
        if "vr_mollweide" in plot_type:        
            pos_galactic = kinematics1.pos_cartesian_to_galactic()
            vel_galactic = kinematics1.vel_cartesian_to_galactic()
            
            pos_galactic2 = kinematics2.pos_cartesian_to_galactic()
            vel_galactic2 = kinematics2.vel_cartesian_to_galactic()

            figname1 = "../../plots/exploration/{}_vr_stars_satellite_faceon_{:03d}.png".format(sim, k)
            fig_title = "{} satellite vr stars; {}-{} kpc; t={:.2f}  Gyr".format(sim, 50, 300, times[k])

            pl.mollweide_projection(pos_galactic[0]*180/np.pi, pos_galactic[1]*180/np.pi, 0, 0, title=fig_title, 
                bmin=-50, bmax=50, nside=40, smooth=1, q=vel_galactic[0], figname=figname1)

            figname1 = "../../plots/exploration/{}_vr_stars_satellite_faceon_300_600_{:03d}.png".format(sim, k)
            fig_title = "{} satellite vr stars; {}-{} kpc; t={:.2f}  Gyr".format(sim, 300, 600, times[k] )
            pl.mollweide_projection(pos_galactic2[0]*180/np.pi, pos_galactic2[1]*180/np.pi, 0, 0, title=fig_title, 
                bmin=-50, bmax=50, nside=40, smooth=1, q=vel_galactic2[0], figname=figname1)

        if  "cartesian_projection" in plot_type:
            print('-> Plotting cartesian projection')
            figname = "../../plots/exploration/{}_DM_stars_projection_300_600_{:03d}.png".format(sim, k)
            pl.multipanel_plot(hfaceon, hfaceon, satellite_faceon, k, sim, figname)

        elif plot_type == 'rho_mollweide':
            print("TBD")

        elif plot_type == "orbital_poles":
            print("TBD")
            
