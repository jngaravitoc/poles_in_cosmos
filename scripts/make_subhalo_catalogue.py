import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from astropy import units as u
import pynbody

sys.path.append("./src/")

sys.path.append("/mnt/home/ecunninagham/python")

import io_gizmo_pynbody as ga
import sys
#import halo_analysis as halo

# 
import pynbody_routines as pr 
import io_gizmo_pynbody as fa
#import plotting as pl

from scipy.linalg import norm
import h5py
import itertools


def write_group(cat_name, snap, subhalos, mask):
    hf = h5py.File(cat_name, 'r+')
    snap_group = hf.create_group('{}'.format(snap))
    snap_group['darkmass'] = subhalos.dark['mass'][mask]
    snap_group['pos']= subhalos.dark['pos'][mask]
    snap_group['vel'] = subhalos.dark['vel'][mask]
    snap_group['stmass'] = subhalos.star['mass']
    snap_group['stpos'] = subhalos.star['pos']
    snap_group['stvel'] = subhalos.star['vel']


    return 0

def sim_angmom(sim, snap, cat_name, rm_sat, rmax=700):
    sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim)
    snap_times = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt".format(sim)
    times = np.loadtxt(snap_times, usecols=3)
    m12 = fa.FIRE(sim, remove_satellite=rm_sat, rm_stellar_sat=True)
    sub_not_sat, sat = m12.subhalos_rotated(snap)
    d = np.sqrt(np.sum(sub_not_sat.dark['pos']**2, axis=1))
    mcut = np.where((np.log10(sub_not_sat.dark['mass']) > 7) & (d<rmax))

    
    write_group(cat_name, snap, sub_not_sat, mcut)


    return 0


if __name__ == "__main__":
       
    sim = sys.argv[1]
    rm_sat = int(sys.argv[2])
    snap_i = 300
    snap_f = 600
    #rm_sat = False

    if bool(rm_sat) == False:
        cat_name = '{}_all_rotated_subhalo_cat.h5py'.format(sim)
        print("Building catalogue without satellite's subhalos")
    elif bool(rm_sat) == True:
        cat_name = '{}_rotated_subhalo_cat.h5py'.format(sim)
        print("Building catalogue with satellite's subhalos")

    hf = h5py.File(cat_name, 'w')

    for k in range(snap_i, snap_f):
        sim_angmom(sim, k, cat_name, bool(rm_sat)) 
