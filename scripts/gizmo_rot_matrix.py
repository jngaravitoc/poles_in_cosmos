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
class FIRE:

  def __init__(self, sim, outpath, figname):
    self.sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim)
    subs_path = '/mnt/home/ecunningham/ceph/latte/{}_res7100/massive_stream/dm_inds.npy'.format(sim)
    self.subs_ids = np.load(subs_path)
    self.outpath = outpath
    self.figname = figname

  def cartessian_density_projections(self, snap):
    subhalos = halo.io.IO.read_catalogs('snapshot', 300, self.sim_directory)
    # Tree
    halt = halo.io.IO.read_tree(simulation_directory=self.sim_directory)

    # Read snapshot
    p = ga.io.Read.read_snapshots(['dark', 'star'], 'snapshot', snap, self.sim_directory, 
                              assign_hosts=True, particle_subsample_factor=1, sort_dark_by_id=True)
   
    # Removing satellite substructure
    npart = len(p['dark'].prop('mass'))
    mask_sub=np.ones(npart, dtype=bool)
    mask_sub[self.subs_ids]=0    

    # Make pynbody halo

    hfaceon = pr.pynbody_halo(p, mask_sub)
    hsideon = pr.pynbody_halo(p, mask_sub)
    pynbody.analysis.angmom.faceon(hfaceon, cen=(0,0,0))
    pynbody.analysis.angmom.sideon(hsideon, cen=(0,0,0))


    #subhalos

    hsub = pr.pynbody_subhalos(subhalos)
    hsub_faceon = pr.pynbody_subhalos(subhalos)

    # Satellite orbit
    sat_id = np.argsort(hsub.dark['mass'])[-2]
    sat_tree_id = subhalos['tree.index'][sat_id]
    satellite = fa.return_tracked_pos(halt, sat_tree_id, pynbody_halo=True)
    satellite_faceon = satellite


    h_rotations = pr.pynbody_halo(p)
    faceon, edgeon = pr.make_pynbody_rotations(h_rotations)

    pynbody.transformation.transform(hsub_faceon, faceon)

    pynbody.transformation.transform(satellite_faceon, faceon)
    figname = self.outpaht + "2d_projection_{:s}".format(snap) + "_.png" 
    pl.multipanel_plot(hfaceon, hsideon, satellite_faceon, snap, figname)




if __name__ == "__main__":
    
    
    
    
    snap_init = 300
    snap_final = 450
    
    sim='m12b'
    figname = "test_projections"
    
    m12b = FIRE(sim, '../plots/exploration/', 'test')
    #m12b.cartessian_density_projections(snap_num, figname)

    
    # Halo catalogue
    
    sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim)
    m12b_subhalos = halo.io.IO.read_catalogs('snapshot', 300, sim_directory)
    # Tree
    halt = halo.io.IO.read_tree(simulation_directory=sim_directory)

    #p0 = ga.io.Read.read_snapshots(['dark', 'star'], 'snapshot', 385, sim_directory, 
    #                              assign_hosts=True, particle_subsample_factor=1)

    subs_path = '/mnt/home/ecunningham/ceph/latte/m12b_res7100/massive_stream/dm_inds.npy'
    subs_ids = np.load(subs_path)                                     
    # load particle data
    
    for k in range(snap_init, snap_final):
        p = ga.io.Read.read_snapshots(['dark', 'star'], 'snapshot', k, sim_directory, 
                                  assign_hosts=True, particle_subsample_factor=1, sort_dark_by_id=True)
       
        # Removing subhalo particles
        # Make pynbody halo

        #subhalos

        hsub = pr.pynbody_subhalos(m12b_subhalos)
        hsub_faceon = pr.pynbody_subhalos(m12b_subhalos)

        # Satellite orbit
        sat_id = np.argsort(hsub.dark['mass'])[-2]
        sat_tree_id = m12b_subhalos['tree.index'][sat_id]
        satellite = fa.return_tracked_pos(halt, sat_tree_id, pynbody_halo=True)
        satellite_faceon = satellite


        h_rotations = pr.pynbody_halo(p)
        faceon, edgeon = pr.make_pynbody_rotations(h_rotations)

        pynbody.transformation.transform(hsub_faceon, faceon)

        pynbody.transformation.transform(satellite_faceon, faceon)

     
