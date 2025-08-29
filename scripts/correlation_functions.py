"""

Author: Nico Garavito-Camargo
Github: jngaravitoc

"""
#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import Angle
import sys
import pynbody

plt.style.use('~/matplotlib.mplstyle')
import gizmo_analysis as ga
import halo_analysis as halo
import nba

import pynbody_routines  as pr 
import FIRE_analysis as fa
import plotting as pl
from gizmo_pynbody_analysis_workflow import FIRE
import analysis as an

def poles_subhalos(snap, rmin=20, rmax=400, satellites=False):
    """
    Compute orbital poles for all the subhalos witin a given range of galactocentric distances. Rmin and Rmax.
    """
    
    f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
    m12_halo = halo.io.IO.read_catalogs('snapshot', snap, sim_directory)
    dist = np.sqrt(np.sum(m12_halo['host.distance']**2, axis=1))
    rcut = np.where((dist>rmin) & (dist<rmax))
                    
    m12_kin = nba.kinematics.Kinematics(m12_halo['host.distance'][rcut], m12_halo['host.velocity'][rcut]*f)
    l, b = m12_kin.orbpole()

    lpol = Angle(l * u.deg)
    lpolw = lpol.wrap_at(360 * u.deg).degree  
    
    if satellites == True :
        stellar_subhalos = m12_halo['star.mass'][rcut]!=-1

        return lpolw[stellar_subhalos], b[stellar_subhalos]
    else :
        return lpolw, b


# Get all the subhalos
def get_halo_satellite(sim, mass_rank):
    sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim)
    m12_subhalos = halo.io.IO.read_catalogs('snapshot', 300, sim_directory)
    halt = halo.io.IO.read_tree(simulation_directory=sim_directory)
    hsub = pr.pynbody_subhalos(m12_subhalos)
    sat_id = np.argsort(hsub.dark['mass'])[mass_rank]
    sat_tree_id = m12_subhalos['tree.index'][sat_id]
    satellite = fa.return_tracked_pos(halt, sat_tree_id, pynbody_halo=True)
    f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
    m12_sat = nba.kinematics.Kinematics(satellite['pos'], satellite['vel']*f)
    l, b = m12_sat.orbpole()
    lpol = Angle(l * u.deg)
    lpolw = lpol.wrap_at(360 * u.deg).degree  
    return lpolw, b


def plot_2dcorrfunc(w, w0, t0, t1, title, figname, hlines=[],  vmin=-0.1, vmax=0.1):
    fig, ax = plt.subplots(1, 1, figsize=(5,4))
    if type(w0) == np.ndarray :
        im = plt.imshow((w+1)/(w0+1) - 1, origin='lower', extent=[0, 180, t0, t1],
                    vmin=vmin, vmax=vmax, aspect='auto', cmap='Spectral')
        cbar = plt.colorbar(im)
        cbar.set_label(r'$\tilde{\omega} (\theta)$')
    else : 
        im = plt.imshow(w, origin='lower', extent=[0, 180, t0, t1],
                    vmin=vmin, vmax=vmax, aspect='auto', cmap='Spectral')
        cbar = plt.colorbar(im)
        cbar.set_label(r'$\omega (\theta)$')
        
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$t\ \rm{[Gyr}]$')
    ax.set_title(title)
    ax.set_xticks([0, 60, 120, 180])
    for n in range(len(hlines)):
        ax.axhline(hlines[n], ls='--', c='k', lw=1)

    plt.savefig(figname, bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    
    
    snap_init = 385
    snap_final = 386
    tsteps = snap_final - snap_init
    nbins = 30
    #nbins = 60
    rmin=300
    rmax=600
    sim='m12c'
    auto = False
    sats = True
    
    wmatrix = np.zeros((tsteps, nbins))
    wmatrix_s = np.zeros((tsteps, nbins))
    
   
    sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim)
    lOP, bOP = poles_subhalos(snap_init, rmin, rmax, satellites=sats)

    snap_times = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt".format(sim)
    times = np.loadtxt(snap_times, usecols=3)
    
    if sim=='m12b':
      lsat, bsat = get_halo_satellite(sim, -2)
      ltimes = [times[385], times[449]]
      
    elif sim=='m12c':
      lsat, bsat = get_halo_satellite(sim, -4)
      ltimes = [times[549]]

    elif sim=='m12f':
      lsat, bsat = get_halo_satellite(sim, -4)
      ltimes = [times[320], times[462]]

    elif sim=='m12i':
      lsat, bsat = get_halo_satellite(sim, -11)
      ltimes = []

    elif sim=='m12m':
      lsat, bsat = get_halo_satellite(sim, -19)
      ltimes = [times[444], times[558]]

    elif sim=='m12r':
      lsat, bsat = get_halo_satellite(sim, -2)
      lsat2, bsat2 = get_halo_satellite(sim, -3)
      lsat3, bsat3 = get_halo_satellite(sim, -5)
      ltimes = [times[477], times[515], times[560]]

    elif sim=='m12w':
      lsat, bsat= get_halo_satellite(sim, -3)
      lsat2, bsat2= get_halo_satellite(sim, -7)
      lsat3, bsat3= get_halo_satellite(sim, -8)
      ltimes = [times[311], times[358], times[490]]



    bins, w0s = an.compute_2d_corrf(lOP, bOP, np.array([np.nanmean(lsat)]), np.array([np.nanmean(bsat)]), nbins)
    bins, w0 = an.compute_2d_corrf(lOP, bOP, np.array([0]), np.array([0]), nbins)
    wmatrix[0] = w0
    wmatrix_s[0] = w0s

    for k in range(snap_init+1, snap_final, 1):
        lOP, bOP = poles_subhalos(k, rmin, rmax, satellites=sats)
        # 300 is the snap @ where the lsat & bsat array starts
        bins, wmatrix_s[k-snap_init] = an.compute_2d_corrf(lOP, bOP, np.array([np.nanmean(lsat)]), 
                                                        np.array([np.nanmean(bsat)]), nbins)
        bins, wmatrix[k-snap_init] = an.compute_2d_corrf(lOP, bOP, np.array([0]), 
                                                        np.array([0]), nbins)
    if auto == True :
      np.savetxt('{}_wmatrix_corrfunc_sat_{}_{}_subhalos_sats_{}.txt'.format(sim, rmin, rmax, sats), wmatrix)
      plot_2dcorrfunc(wmatrix_s, 0, times[snap_init], times[snap_final],  r'${}$'.format(sim), '{}_2d_corrfunc_sat_{}_{}_{}.pdf'.format(sim, sats, rmin, rmax), ltimes, vmin=-2, vmax=2)
    else :
      np.savetxt('{}_wmatrix_corrfunc_{}_{}_subhalos_sats_{}.txt'.format(sim, rmin, rmax, sats), wmatrix)
      plot_2dcorrfunc(wmatrix, wmatrix[0], times[snap_init], times[snap_final],  r'${}$'.format(sim), '{}_2d_corrfunc_{}_{}_{}.pdf'.format(sim, sats, rmin, rmax), ltimes, vmin=-0.1, vmax=0.1)
