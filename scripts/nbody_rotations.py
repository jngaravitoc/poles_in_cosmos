"""
Routines for defining pynbody halos and perfoming rotations 

Dependencies:
  - pynbody

author: Nico Garavito-Camargo
github: jngaravitoc

"""


import numpy as np
import pynbody
from pynbody import filt

def pynbody_halo(particles):
    """
    TODO:
    - Check units 
    - Do we need gas? if so it can be added as:
        halo_pynb.gas['pos'] = particles['gas_pos']
    """

    ndark = len(particles['dm_mass'])
    nstar = len(particles['star_mass'])
    halo_pynb = pynbody.new(dark=int(ndark), star=int(nstar), order='dark,star')
    
    halo_pynb.dark['pos'] = particles['dm_pos']
    halo_pynb.dark['vel'] = particles['dm_vel']
    halo_pynb.dark['mass'] = particles['dm_mass']

    halo_pynb.star['pos'] = particles['star_pos']
    halo_pynb.star['vel'] = particles['star_vel']
    halo_pynb.star['mass'] = particles['star_mass']

    halo_pynb.dark['pos'].units = 'kpc'
    halo_pynb.dark['vel'].units = 'km s**-1'
    halo_pynb.dark['mass'].units = '1e10 Msol'

    halo_pynb.star['pos'].units = 'kpc'
    halo_pynb.star['vel'].units = 'km s**-1'
    halo_pynb.star['mass'].units = '1e10 Msol'
    
    return halo_pynb

def make_pynbody_rotations(halo):
    cen = halo[pynbody.filt.Sphere("5 kpc")]
    Lh = pynbody.analysis.angmom.ang_mom_vec(cen)
    Tx_faceon = pynbody.analysis.angmom.calc_faceon_matrix(Lh)
    Tx_sideon = pynbody.analysis.angmom.calc_sideon_matrix(Lh)
    return Tx_faceon, Tx_sideon



if __name__ == "__main__":
    # Read data
    dm_pos = np.zeros((10,3))
    dm_pos[:,0] = np.arange(0, 10)
    dm_vel = dm_pos
    dm_mass = np.ones(10)
    star_pos = np.ones((10,3))
    star_vel = np.ones((10,3))
    star_mass = np.ones(10)


    # Make dictionary for pynbody
    sim_data = { "dm_pos": dm_pos,
                 "dm_vel": dm_vel,
                 "dm_mass": dm_mass,
                 "star_pos": star_pos,
                 "star_vel": star_vel,
                 "star_mass": star_mass,
    }
    # Make pynbody halo
    pynb_halo = pynbody_halo(sim_data)
    
    faceon, edgeon = make_pynbody_rotations(pynb_halo)
    # Rotate 
    PA_rotation = pynbody.transformation.transform(pynb_halo, faceon)
    # Add LMC

    # Reverse rotation 
    PA_rotation.revert()

