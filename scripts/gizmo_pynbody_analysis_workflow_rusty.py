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
    stars_path = '/mnt/home/ecunningham/ceph/latte/{}_res7100/massive_stream/new_z0_inds.npy'.format(sim)
    self.stats_ids = np.load(stars_path)
    self.outpath = outpath
    self.figname = figname
    self.sim = sim

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
    figname = self.outpath + "2d_projection_{}".format(snap) + "_.png" 
    pl.multipanel_plot(hfaceon, hsideon, satellite_faceon, snap, self.sim,  self.outpath + self.figname)




if __name__ == "__main__":
    
    
    
    snap = int(sys.argv[1])    
    snap_base = 300
    assert snap >= snap_base

    
    sim='m12c'
    figname = '../plots/exploration/m12b_DM_stars_projections'
    
    m12b = FIRE(sim, '../plots/exploration/', '{}_DM_stars_projections'.format(sim))
     
    
    # Halo catalogue
    
    sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim)
    m12b_subhalos = halo.io.IO.read_catalogs('snapshot', snap_base, sim_directory)
    # Tree
    halt = halo.io.IO.read_tree(simulation_directory=sim_directory)


    subs_path = '/mnt/home/ecunningham/ceph/latte/{}_res7100/massive_stream/dm_inds.npy'.format(sim)
    subs_ids = np.load(subs_path)                                     
    # load particle data
    
    #for k in range(snap_init, snap_final):
    
    m12b.cartessian_density_projections(snap)
    p = ga.io.Read.read_snapshots(['dark', 'star'], 'snapshot', snap, sim_directory, 
                              assign_hosts=True, particle_subsample_factor=1, sort_dark_by_id=True)
   
    # Removing subhalo particles
    
    npart = len(p['dark'].prop('mass'))
    
    # remove subhalo 
    mask_sub=np.ones(npart, dtype=bool)
    mask_sub[subs_ids]=0 
    
    # Only leave subhalo
    #mask_sub=np.zeros(npart, dtype=bool)
    #mask_sub[subs_ids]=1 
    
    # Make pynbody halo
    hfaceon = pr.pynbody_halo(p, mask_sub)
    hsideon = pr.pynbody_halo(p, mask_sub)
    pynbody.analysis.angmom.faceon(hfaceon, cen=(0,0,0))
    pynbody.analysis.angmom.sideon(hsideon, cen=(0,0,0))


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

    #multipanel_plot(hfaceon, hsideon, satellite_faceon, snap, figname)
    
    pos_dm = hfaceon.dark['pos']
    f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
    vel_dm = hfaceon.dark['vel']*f
    dist_dm = np.sqrt(np.sum(pos_dm**2, axis=1))
    dist_cut = np.where((dist_dm> 50) & (dist_dm< 300)) 
    dm_kinematics = nba.kinematics.Kinematics(pos_dm[dist_cut],  vel_dm[dist_cut])
    dm_l_host, dm_b_host = dm_kinematics.pos_cartesian_to_galactic()
    dm_OP_l_host, dm_OP_b_host = dm_kinematics.orbpole()
    

    sat_kinematics = nba.kinematics.Kinematics(satellite_faceon.dark['pos'],  satellite_faceon.dark['vel'])
    sat_l_host, sat_b_host = sat_kinematics.pos_cartesian_to_galactic()
    sat_OP_l_host, sat_OP_b_host = sat_kinematics.orbpole()

    dm_figname = "../plots/exploration/{}_galacitc_faceon_{:03d}.png".format(sim, snap)
    times = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt'.format(sim)
    t_snap = np.loadtxt(times, usecols=3)
    #title_dm = "{} DM; {}-{} kpc; t={:.2f}  Gyr".format('m12b', 50, 300,  t_snap[k] )
    #pl.mollweide_projection(dm_l_host*180/np.pi, dm_b_host*180/np.pi, sat_l_host[k-300]*180/np.pi, sat_b_host[k-300]*180/np.pi, 
    #                     title=title_dm, bmin=100, bmax=1000,
    #                     nside=40, figname=dm_figname)
                         
    title_dm = "{} satellite poles DM; {}-{} kpc; t={:.2f}  Gyr".format(sim, 50, 300, t_snap[k] )
    dm_figname =  "../plots/exploration/{}_OP_satellite_faceon_{:03d}.png".format(sim, snap)
    pl.mollweide_projection(dm_OP_l_host, dm_OP_b_host, sat_OP_l_host[snap-snap_base], sat_OP_b_host[snap-snap_base], 
                         title=title_dm, bmin=100, bmax=1000,
                         nside=40, smooth=1, figname=dm_figname)
 
