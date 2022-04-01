"""
Script to compute orbital poles maps in m12b Fire halo
author: ecunningham, jngaravitoc
02/2022 -

TODO:
-----
- Plot orbital pole of the satellite
- In which direction the COM off-set should point too in the orbital poles map.

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
import nba
import healpy as hp
from  healpy.newvisufunc import projview, newprojplot
from scipy.linalg import norm


def get_3d_rotation(vec1, vec2):
    """
    Computes the rotation matrix from going to vec1 to vec2
    
    """
    assert np.dot(vec1, vec2) != -1, "Method not valid for vec1 = -vec2"
    v1xv2 = np.cross(vec1, vec2)
    v1dv2 = np.dot(vec1, vec2)
    s = norm(v1xv2)
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    V = np.array([[0, -v1xv2[2], v1xv2[1]], [v1xv2[2], 0, -v1xv2[0]], [-v1xv2[1], v1xv2[0], 0]])
    R = I + V + (np.dot(V, V) / (1+v1dv2))
    return R

def cartessian_projection(pos, orbit, figname):
    fig, ax = plt.subplots(1, 2, figsize=(10,4))
    ax[0].hist2d(pos[:,0], pos[:,1],  bins=np.linspace(-100,100,800), norm=LogNorm())
    ax[1].hist2d(pos[:,1], pos[:,2],  bins=np.linspace(-100,100,800), norm=LogNorm())
    ax[0].set_xlabel("x[kpc]")
    ax[0].set_ylabel("y[kpc]")
    ax[1].set_xlabel("y[kpc]")
    ax[1].set_ylabel("z[kpc]")
    ax[0].plot(orbit[:,0], orbit[:,1], c='r')
    ax[1].plot(orbit[:,1], orbit[:,2], c='r')
    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    return 0

def mollweide_projection(l, b, l2, b2, rmin, rmax, bmin, bmax, nside, figname, smooth=5):

    """
    Makes mollweide plot using healpix
    Parameters:
    ----------- 
    l : numpy.array
    b : numpy.array
    """
 
    mwlmc_indices = hp.ang2pix(nside,  (90-b)*np.pi/180., l*np.pi/180.)
    npix = hp.nside2npix(nside)
 
    idx, counts = np.unique(mwlmc_indices, return_counts=True)
    degsq = hp.nside2pixarea(nside, degrees=True)
    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts/degsq
    map_smooth = hp.smoothing(hpx_map, fwhm=smooth*np.pi/180)
    
  
    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    #twd_map = hp.mollview(np.array([map_smooth]), rot=(180, 90, 0), return_projected_map=True, min=0, max=500)
    plt.close()
    projview(
    hpx_map,
    coord=["G"],
    graticule=True,
    graticule_labels=True,
    rot=(120, 0, 0),
    unit=" ",
    xlabel="Galactic Longitude (l) ",
    ylabel="Galactic Latitude (b)",
    cb_orientation="horizontal",
    min=bmin,
    max=bmax,
    latitude_grid_spacing=45,
    projection_type="mollweide",
    title=" {}-{} kpc".format(int(rmin),int(rmax)),)
	
    #newprojplot(theta=np.radians(90-(b2)), phi=np.radians(l2-120), marker="*", color="r", markersize=5 )
    newprojplot(theta=np.radians(90-(b2[0])), phi=np.radians(l2[0]-120), marker="*", color="r", markersize=5 )
    #newprojplot(theta=np.radians(90-(b2[1])), phi=np.radians(l2[1]-120), marker="*", color="w", markersize=2 )

    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    return 0




	
def rotatate_halo(pos, vel):

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
    print("rotating halo with angles: ", x_ang, y_ang, z_ang)
    part_rot=ut.coordinate.get_coordinates_rotated(pos, rotation_angles=[x_ang, 0., z_ang])
    halo_pos_rot=ut.coordinate.get_coordinates_rotated(halo_pos, rotation_angles=[x_ang, 0., z_ang])

   
    vel_part_rot=ut.coordinate.get_coordinates_rotated(vel, rotation_angles=[x_ang, 0., z_ang])
    halo_pos_rot=ut.coordinate.get_coordinates_rotated(halo_pos, rotation_angles=[x_ang, 0., z_ang])


    # Satellite orbital poles
    sat_vel=np.array([-181.44563,262.926,141.28804])*to_kpc_gyr
    sat_pos=np.array([32.908375,16.943731,8.137191])
    sat_pos_rot=ut.coordinate.get_coordinates_rotated(sat_pos, rotation_angles=[x_ang, 0., z_ang])
    sat_vel_rot=ut.coordinate.get_coordinates_rotated(sat_vel, rotation_angles=[x_ang, 0., z_ang])
    print("Satellite position and velocity", sat_pos, sat_vel)


    return part_rot, vel_part_rot, halo_pos_rot, sat_pos_rot, sat_vel_rot

## Orbital poles routines



if __name__ == "__main__":
    # Parameters

    simulation_directory = '/mnt/ceph/users/firesims/fire2/metaldiff/m12b_res7100'
    snap_number = 385
    outpath = "./"
   

    part = ga.io.Read.read_snapshots(['star', 'dark'],'snapshot', snap_number, simulation_directory,
        assign_hosts=True, assign_orbits=True, assign_hosts_rotation=True)

    to_kpc_gyr = ((1*u.km/u.s).to(u.kpc/u.Gyr)).value
	# * Stars
    stars_host_distance=part['star'].prop('host.distance.principal')
    stars_host_velocity=part['star'].prop('host.velocity.principal')*to_kpc_gyr
    not_in_subs_stars=np.loadtxt('/mnt/ceph/users/ecunningham/latte/m12b_{:0>3d}_unbound_star_indices.txt'.format(snap_number))
    in_subs_stars=np.ones(len(part['star']['potential']), dtype='bool')
    in_subs_stars[not_in_subs_stars.astype(int)]=False
    # Density plot of dark matter particles withoug substrucure

    # *  Rotate halo without substructure
    #pos_part_rot, vel_part_rot, halo_pos_rot, sat_pos, sat_vel = rotatate_halo(dm_host_distance[not_in_subs.astype(int)], dm_host_velocity[not_in_subs.astype(int)])
    # Rotate halo with substructure 
    #pos_part_rot, vel_part_rot, halo_pos_rot, sat_pos, sat_vel = rotatate_halo(dm_host_distance, dm_host_velocity)
    #print("N total particles, N particles w/o subs", len(pos_part_rot), len(dm_host_distance))

    # *  Dark matter 
    not_in_subs_dark=np.load('/mnt/ceph/users/ecunningham/latte/m12b_{:0>3d}_unbound_dark_indices.npy'.format(snap_number))
    in_subs_dark=np.ones(len(part['dark']['potential']), dtype='bool')
    in_subs_dark[not_in_subs_dark.astype(int)]=False
    
    dm_host_distance=part['dark'].prop('host.distance.principal')
    dm_host_velocity=part['dark'].prop('host.velocity.principal')*to_kpc_gyr
    # Satellites orbit:
    #future_orb=np.load('/mnt/home/ecunningham/ceph/latte/m12b_385_lmc_future_orbit.npy')
    previous_orb=np.load('/mnt/home/ecunningham/ceph/latte/m12b_385_lmc_prev_orbit.npy')
    R = part.host['rotation']
    #future_orb_rot=np.dot(R[0], future_orb.T).T
    previous_orb_rot=np.dot(R[0], previous_orb.T).T
    #cartessian_projection(dm_host_distance[not_in_subs_dark.astype(int)], np.zeros((3,3)), "rho2_dark_particles_cartessian.png")
    #cartessian_projection(dm_host_distance[in_subs_dark], np.zeros((3,3)), "rho2_dark_subs_particles_cartessian.png")
    #artessian_projection(stars_host_distance, previous_orb_rot, "rho_stellar_particles_cartessian_320.png")


    sat_vel=np.array([-181.44563,262.926,141.28804])*to_kpc_gyr
    sat_pos=np.array([32.908375,16.943731,8.137191])
    sat_pos_rot = np.dot(R[0], sat_pos.T).T
    sat_vel_rot = np.dot(R[0], sat_vel.T).T
	#mass =  part['star']['mass'][not_in_subs]
    # Re-center halo
    #rcom, vcom = nba.com.shrinking_sphere(pos_part_rot, vel_part_rot, mass)
    #print("New rcom", rcom)
    #print("New vcom", vcom)


    #pos_com = nba.com.re_center(pos_part_rot, rcom)
    #vel_com = nba.com.re_center(vel_part_rot, vcom)

    sat_kinematics = nba.kinematics.Kinematics(np.array([sat_pos_rot]), np.array([sat_vel_rot]))
    op_l_sat, op_b_sat = sat_kinematics.orbpole()
    l_sat, b_sat = sat_kinematics.pos_cartesian_to_galactic()
    print("Satellite orb poles", op_l_sat, op_b_sat)
   
    orbit_kinematics = nba.kinematics.Kinematics(previous_orb_rot, previous_orb_rot)
    l_orbit, b_orbit = orbit_kinematics.pos_cartesian_to_galactic()
    
    
    dm_dist = np.sqrt(np.sum(dm_host_distance[not_in_subs_dark.astype(int)]**2, axis=1))
    dm_dist_subs = np.sqrt(np.sum(dm_host_distance[in_subs_dark]**2, axis=1))
    print("len of subsparticles", len(dm_dist_subs))
    #dm_dist = np.sqrt(np.sum(dm_host_distance**2, axis=1))
    stars_dist = np.sqrt(np.sum(stars_host_distance[not_in_subs_stars.astype(int)]**2, axis=1))
    #stars_dist = np.sqrt(np.sum(stars_host_distance**2, axis=1))
    r = np.arange(0, 551, 50)
    #r = np.arange(0, 51, 10)

    for k in range(len(r)-1):
    	dm_dist_cut = np.where((dm_dist<r[k+1]) & (dm_dist>r[k]))
    	#dm_dist_cut_subs = np.where((dm_dist_subs<r[k+1]) & (dm_dist_subs>r[k]))
    	#stars_dist_cut = np.where((stars_dist<r[k+1]) & (stars_dist>r[k]))
    	#dm_m12b_kinematics_2 = nba.kinematics.Kinematics(dm_host_distance[not_in_subs_dark.astype(int)][dm_dist_cut], dm_host_velocity[not_in_subs_dark.astype(int)][dm_dist_cut])
    	#dm_m12b_kinematics_subs = nba.kinematics.Kinematics(dm_host_distance[in_subs_dark][dm_dist_cut_subs], dm_host_velocity[in_subs_dark][dm_dist_cut_subs])
    	dm_m12b_kinematics_2 = nba.kinematics.Kinematics(dm_host_distance[dm_dist_cut], dm_host_velocity[dm_dist_cut])
    	#stars_m12b_kinematics_2 = nba.kinematics.Kinematics(stars_host_distance[not_in_subs_stars.astype(int)][stars_dist_cut], stars_host_velocity[not_in_subs_stars.astype(int)][stars_dist_cut])
    	#dm_op_l_2, dm_op_b_2 = dm_m12b_kinematics_2.orbpole()
    	#dm_op_l_subs, dm_op_b_subs = dm_m12b_kinematics_subs.orbpole()
    	dm_l_host, dm_b_host = dm_m12b_kinematics_2.pos_cartesian_to_galactic()
    	#stars_op_l_2, stars_op_b_2 = stars_m12b_kinematics_2.orbpole()
    	#stars_l_host, stars_b_host = stars_m12b_kinematics_2.pos_cartesian_to_galactic()
    	dm_figname = "orb_poles_m12b_dark_particles_shell_{:03d}.png".format(k)
    	dm_figname_subs = "orb_poles_m12b_dark_particles_subs_shell_{:03d}.png".format(k)
    	#stars_figname = "orb_poles_m12b_star_particles_shell_{:03d}.png".format(k)
    	#mollweide_projection(dm_op_l_2, dm_op_b_2, [op_l_sat], [op_b_sat], r[k], r[k+1], bmin=0, bmax=700, nside=60, figname=dm_figname)
    	#mollweide_projection(dm_op_l_subs, dm_op_b_subs, [op_l_sat], [op_b_sat], r[k], r[k+1], bmin=0, bmax=700, nside=60, figname=dm_figname_subs)
    	#mollweide_projection(stars_op_l_2, stars_op_b_2, op_l_sat, op_b_sat, r[k], r[k+1], bmin=0, bmax=10, nside=30, figname=stars_figname)
    	dm_figname = "rho_m12b_dark_particles_shell_{:03d}.png".format(k)
    	#stars_figname = "rho_m12b_star_particles_shell_{:03d}.png".format(k)
    	mollweide_projection(dm_l_host*180/np.pi, dm_b_host*180/np.pi, [l_orbit*180/np.pi], [b_orbit*180/np.pi], r[k], r[k+1], bmin=0, bmax=700, nside=40, figname=dm_figname)
    	#mollweide_projection(stars_l_host*180/np.pi, stars_b_host*180/np.pi, 0*180/np.pi, 0*180/np.pi, r[k], r[k+1],bmin=0, bmax=7, nside=40, figname=stars_figname)

    
		
