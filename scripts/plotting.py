"""
Plotting routines used for the orbital poles project

Dependencies:
  - Matplitlib
  - Healpy
  - pynbody 

author: Nico Garavito-Camargo
github: jngaravitoc



"""


import numpy as np
import matplotlib.pyplot as plt
import pynbody
import sys
from pynbody import filt
from astropy import units as u
import nba
import healpy as hp
from healpy.newvisufunc import projview, newprojplot

sys.path.append("/mnt/home/ecunningham/python")
plt.style.use('~/matplotlib.mplstyle')

    
def multipanel_plot(hf, satellite_faceon, snap, sim, figname):
    times = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt'.format(sim)
    t_snap = np.loadtxt(times, usecols=3)

   

    fig, ax = plt.subplots(2, 3, figsize=(20, 13))
    dmproj = pynbody.plot.image(hf.d, width=1200, cmap='Greys', subplot=ax[0][0], show_cbar=False, av_z=True)
    dmproj2 = pynbody.plot.image(hs.d, width=1200, cmap='Greys', subplot=ax[1][0], show_cbar=False, av_z=True)

    dmproj = pynbody.plot.image(hf.s, width=300, cmap='inferno', subplot=ax[0][1], show_cbar=False, av_z=True)
    dmproj2 = pynbody.plot.image(hs.s, width=300, cmap='inferno', subplot=ax[1][1], show_cbar=False, av_z=True)

    dmproj = pynbody.plot.image(hf.s, width=100, cmap='inferno', subplot=ax[0][2], show_cbar=False, av_z=True)
    dmproj2 = pynbody.plot.image(hs.s, width=100, cmap='inferno', subplot=ax[1][2], show_cbar=False, av_z=True)


    ax[0][0].plot(satellite_faceon.dark['pos'][:snap-300,0], satellite_faceon.dark['pos'][:snap-300,1], c='k', ls='--', alpha=0.6, lw=1)
    ax[1][0].plot(satellite_faceon.dark['pos'][:snap-300,0], satellite_faceon.dark['pos'][:snap-300,2], c='k', ls='--', alpha=0.6, lw=1)

    ax[0][1].plot(satellite_faceon.dark['pos'][:snap-300,0], satellite_faceon.dark['pos'][:snap-300,1], c='w', ls='--', alpha=0.6, lw=1)
    ax[1][1].plot(satellite_faceon.dark['pos'][:snap-300,0], satellite_faceon.dark['pos'][:snap-300,2], c='w', ls='--', alpha=0.6, lw=1)

    ax[0][2].plot(satellite_faceon.dark['pos'][:snap-300,0], satellite_faceon.dark['pos'][:snap-300,1], c='w', ls='--', alpha=0.6, lw=1)
    ax[1][2].plot(satellite_faceon.dark['pos'][:snap-300,0], satellite_faceon.dark['pos'][:snap-300,2], c='w', ls='--', alpha=0.6, lw=1)



    ax[0][1].set_ylabel('')
    ax[0][2].set_ylabel('')
    ax[1][1].set_ylabel('')
    ax[1][2].set_ylabel('')

    ax[0][0].set_xlabel('$x\mathrm{[kpc]}$')
    ax[0][1].set_xlabel('$x\mathrm{[kpc]}$')
    ax[0][2].set_xlabel('$x\mathrm{[kpc]}$')

    ax[1][0].set_xlabel('$x\mathrm{[kpc]}$')
    ax[1][1].set_xlabel('$x\mathrm{[kpc]}$')
    ax[1][2].set_xlabel('$x\mathrm{[kpc]}$')


    ax[0][0].set_ylabel('$y\mathrm{[kpc]}$')
    ax[1][0].set_ylabel('$z\mathrm{[kpc]}$')

    ax[0][0].set_title('$\mathrm{Dark\ Matter}$', fontsize=16)
    ax[0][1].set_title('$\mathrm{Stars\ outer\ halo}$', fontsize=16)
    ax[0][2].set_title('$\mathrm{Stars\ inner\ halo}$', fontsize=16)
    fig.suptitle('$t={:.2f}$'.format(t_snap[snap]) + r' $\rm{Gyr;}\ \rm{Sim:\ }$ ' + '{}'.format(sim), fontsize=18, y=0.95)

    ax[0][0].text(-500, 400, r"$\rm{Face-on}$", c='k', fontsize=18)
    ax[1][0].text(-500, 400, r"$\rm{Edge-on}$", c='k', fontsize=18)

    ax[0][0].set_xlim(-600, 600)
    ax[1][0].set_xlim(-600, 600)

    ax[0][1].set_xlim(-150, 150)
    ax[1][1].set_xlim(-150, 150)

    ax[0][2].set_xlim(-50, 50)
    ax[1][2].set_xlim(-50, 50)

    ax[0][0].set_ylim(-600, 600)
    ax[1][0].set_ylim(-600, 600)

    ax[0][1].set_ylim(-150, 150)
    ax[1][1].set_ylim(-150, 150)

    ax[0][2].set_ylim(-50, 50)
    ax[1][2].set_ylim(-50, 50)
    plt.savefig(figname + "_{:03d}.png".format(snap), bbox_inches='tight')
    plt.close()
    
    
def mollweide_projection(l, b, l2, b2, title, bmin, bmax, nside, smooth, **kwargs):

    """
    Makes mollweide plot using healpix
    Parameters:
    ----------- 
    l : numpy.array
    b : numpy.array
    """
 
    times = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt'.format('m12b')

    mwlmc_indices = hp.ang2pix(nside,  (90-b)*np.pi/180., l*np.pi/180.)
    npix = hp.nside2npix(nside)
 
    idx, counts = np.unique(mwlmc_indices, return_counts=True)
    degsq = hp.nside2pixarea(nside, degrees=True)
    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts/degsq
    map_smooth = hp.smoothing(hpx_map, fwhm=smooth*np.pi/180)
    
  
    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    plt.close()
    projview(
      map_smooth,
      coord=["G"],
      graticule=True,
      graticule_labels=True,
      rot=(0, 0, 0),
      unit=" ",
      xlabel="Galactic Longitude (l) ",
      ylabel="Galactic Latitude (b)",
      cb_orientation="horizontal",
      min=bmin,
      max=bmax,
      latitude_grid_spacing=45,
      projection_type="mollweide",
      title=title,)
	
    newprojplot(theta=np.radians(90-(b2)), phi=np.radians(l2), marker="*", color="r", markersize=5)
    #newprojplot(theta=np.radians(90-(b2[0])), phi=np.radians(l2[0]-120), marker="*", color="r", markersize=5 )
    #newprojplot(theta=np.radians(90-(b2[1])), phi=np.radians(l2[1]-120), marker="*", color="w", markersize=2 )
    
    if 'figname' in kwargs.keys():
        print("* Saving figure in ", kwargs['figname'])
        plt.savefig(kwargs['figname'], bbox_inches='tight')
        plt.close()
    #return 0

        
