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


"""

#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import sys
import pynbody

import gizmo_analysis as ga
import halo_analysis as halo
import nba

# 
import pynbody_routines  as pr 
from io_gizmo_pynbody  import FIRE
import analysis as an

if __name__ == "__main__":
    
    
    snap_init = 530
    snap_final = 600
    sim='m12b'
    sats = False
    rmin = 50
    rmax = 300
    
    snap_times = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt'.format(sim)
    times = np.loadtxt(snap_times, usecols=3)
    coms_r = np.zeros((snap_final-snap_init, 3))
    coms_v = np.zeros((snap_final-snap_init, 3))
    j=0 
    m12 = FIRE(sim, remove_satellite=True)
    for k in range(snap_init, snap_final, 10):
        # face on particle data halo
        hfaceon, _ = m12.rotated_halo(k)

        
        pos_dm = hfaceon.dark['pos']
        f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
        vel_dm = hfaceon.dark['vel']*f
        dist_dm = np.sqrt(np.sum(pos_dm**2, axis=1))

        dist_cut1 = np.where((dist_dm> rmin) & (dist_dm< rmax)) 
        npart = len(dist_cut1[0])      
        comx, comv = nba.com.com_methods.mean_pos(pos_dm[dist_cut1],  vel_dm[dist_cut1], np.ones(npart))
        print(comx)
        print(comv)
        coms_r[j] = comx
        coms_v[j] = comv
        j+=1
    np.savetxt("test_comr2.txt", coms_r)
    np.savetxt("test_comv2.txt", coms_v)
