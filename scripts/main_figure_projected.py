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
import FIRE_analysis as fa
import plotting as pl


from gizmo_pynbody_analysis_workflow import FIRE

if __name__ == "__main__":
    
    
    snap_init = 450
    snap_final = 451
    
    sim='m12w'
    

    
    sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim)
    
    m12 = FIRE(sim, remove_satellite=False)
    sat_ids = [-3, -7, -8]
    for k in [340, 358, 364]:
    #for k in range(snap_init, snap_final):
        # face on particle data halo
        hfaceon = m12.rotated_halo(k)

        # Satellite orbit
        subhalos_faceon, satellite_faceon = m12.subhalos_rotated(k, sat_index=sat_ids)
      
        
        figname = "../plots/exploration/{}_DM_stars_projection_{:03d}.png".format(sim, k)
        pl.multipanel_plot(hfaceon,  satellite_faceon,  k, sim, figname)

