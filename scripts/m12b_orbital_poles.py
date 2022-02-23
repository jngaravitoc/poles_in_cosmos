"""
Script to compute orbital poles maps in m12b Fire halo
author: ecunningham, jngaravitoc
02/2022 -

"""

#!/usr/bin/env python
# coding: utf-8



import numpy as np
from numpy import Inf

import sys
import schwimmbad
from scipy.optimize import curve_fit
from astropy import units as u
import gala.potential as gp
from gala.units import galactic, solarsystem, dimensionless
from gala.potential import scf

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sklearn as sk

import sys
sys.path.append("/mnt/home/ecunningham/python")
sys.path.append("/mnt/home/nico/codes/nbody_analysis/src")

import gizmo_analysis as ga
import halo_analysis as ra
import utilities as ut


def rotatate_halo(pos):

    #
    halo_pos=np.array([32.908375,16.943731,8.137191])
    halo_r=(np.sum(halo_pos**2)**0.5)
    halo_rp=(np.sum(halo_pos[:2]**2)**0.5)
    z_ang=np.arcsin(halo_pos[2]/halo_r)

    halo_eigen_vec=halo_pos/halo_r
    #print(np.arc)
    x_ang=np.arccos(halo_pos[0]/halo_rp)
    y_ang=np.arcsin(halo_pos[1]/halo_rp)
    ##Note here that x_ang=y_ang=theta; rotation in x-y plane. Rotate in plane and then rotate Satellite down.
    #print(np.arccos(halo_pos[0]/halo_rp))

    part_rot=ut.coordinate.get_coordinates_rotated(pos, rotation_angles=[x_ang, 0., z_ang])
    halo_pos_rot=ut.coordinate.get_coordinates_rotated(halo_pos, rotation_angles=[x_ang, 0., z_ang])

    vel_part_rot=ut.coordinate.get_coordinates_rotated(vel, rotation_angles=[x_ang, 0., z_ang])
    #halo_pos_rot=ut.coordinate.get_coordinates_rotated(halo_pos, rotation_angles=[x_ang, 0., z_ang])


    return part_rot, vel_part, rot, halo_pos_rot

## Orbital poles routines



if __name__ == "__main__":
    # Parameters

    simulation_directory = '/mnt/ceph/users/firesims/fire2/metaldiff/m12b_res7100'
    snap_number = 385
    outpath = "./"
    




    part = ga.io.Read.read_snapshots(['star', 'dark'],'snapshot', snap_number, simulation_directory,
        assign_hosts=True, assign_orbits=True)

    not_in_subs=np.load('/mnt/ceph/users/ecunningham/latte/m12b_{:0>3d}_unbound_dark_indices.npy'.format(snap_number))

    sub_mask=np.zeros(len(part['dark']['potential']), dtype='bool')
    sub_mask[not_in_subs.astype(int)]=True
    # Density plot of dark matter particles withoug substrucure
    dm_host_distance=part['dark'].prop('host.distance')
    dm_host_velocity=part['dark'].prop('host.velocity')

    # Rotate halo
    pos_part_rot, vel_part_rot, halo_pos_rot = rotatate_halo(dm_host_distance[not_in_subs], dm_host_velocity[not_in_subs])
