"""
Scripts to analyze FIRE simulatios
author: jngaravitoc, Emily Cunningham

02/2022 -

TODO:
-----
- Implement parallelization
- Implement arg parser
- Add routine to plot densities

"""

#!/usr/bin/env python
# coding: utf-8



import numpy as np
from numpy import Inf

import sys
import schwimmbad
from scipy.optimize import curve_fit
from astropy import units as u
#import gala.potential as gp
#from gala.units import galactic, solarsystem, dimensionless
#from gala.potential import scf

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#import sklearn as sk

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


## Rotations 

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


## Plotting routines 

def cartessian_projection(pos, figname):
    fig, ax = plt.subplots(1, 2, figsize=(10,4))
    ax[0].hist2d(pos[:,0], pos[:,1],  bins=np.linspace(-100,100,800), norm=LogNorm())
    ax[1].hist2d(pos[:,0], pos[:,2],  bins=np.linspace(-100,100,800), norm=LogNorm())
    ax[0].set_xlabel("x[kpc]")
    ax[0].set_ylabel("y[kpc]")
    ax[1].set_xlabel("x[kpc]")
    ax[1].set_ylabel("z[kpc]")
    #ax[0].plot(orbit[:,0], orbit[:,1], c='r')
    #ax[1].plot(orbit[:,1], orbit[:,2], c='r')
    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    return 0

def mollweide_projection(l, b, l2, b2, rmin, rmax, bmin, bmax, nside,  smooth=5, **kwargs):

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
    #newprojplot(theta=np.radians(90-(b2[0])), phi=np.radians(l2[0]-120), marker="*", color="r", markersize=5 )
    #newprojplot(theta=np.radians(90-(b2[1])), phi=np.radians(l2[1]-120), marker="*", color="w", markersize=2 )
    
    if 'figname' in kwargs.keys():
        print("* Saving figure in ", kwargs['figname'])
        plt.savefig(kwargs['figname'], bbox_inches='tight')
        plt.close()
    return 0



def rotate_halo(pos, vel):

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
    
    pos_part_rot=ut.coordinate.get_coordinates_rotated(pos, rotation_angles=[x_ang, 0., z_ang])
    halo_pos_rot=ut.coordinate.get_coordinates_rotated(halo_pos, rotation_angles=[x_ang, 0., z_ang])

   
    vel_part_rot=ut.coordinate.get_coordinates_rotated(vel, rotation_angles=[x_ang, 0., z_ang])
    halo_pos_rot=ut.coordinate.get_coordinates_rotated(halo_pos, rotation_angles=[x_ang, 0., z_ang])


    # Satellite orbital poles
    to_kpc_gyr = ((1*u.km/u.s).to(u.kpc/u.Gyr)).value
    sat_vel=np.array([-181.44563,262.926,141.28804])*to_kpc_gyr
    sat_pos=np.array([32.908375,16.943731,8.137191])
    sat_pos_rot=ut.coordinate.get_coordinates_rotated(sat_pos, rotation_angles=[x_ang, 0., z_ang])
    sat_vel_rot=ut.coordinate.get_coordinates_rotated(sat_vel, rotation_angles=[x_ang, 0., z_ang])
    print("-> Satellite position and velocity", sat_pos, sat_vel)


    return pos_part_rot, vel_part_rot # halo_pos_rot, sat_pos_rot, sat_vel_rot


## Orbital poles routines

class FIRE_analysis(object):
    def __init__(self, sim, rmin, rmax, part_type, analysis_type,  outpath, figure_name, **kwargs):
        """
            sim: str
                currently only tested with: m21b, m12c, m12f, m12i, m12w for the
                7100 res.
            rmin : float
                Minimum radius at which to slice the halo in kpc
            rmax : float
                Maximum radius at which to slice the halo in kpc
            part_type : list of strings
                Particule type available: 'star', 'dark' if both are desired
                please pass a list of strings ['star', 'dark']
            snap : int
                Snapshot number
            analysis_type : str
                Type of analysis used. 
                Currently supported are: (orbital_poles, density_plots)
            figure_name : str
                Name of the figure
            output_path : str
                Path to where the figure will be save


            **kwargs:

                remove_subs : bool
                    Remove substructure using IDs computed by Emily Cunnigham.



        """
        self.sim = sim
        self.rmin = rmin
        self.rmax = rmax
        self.part_type = part_type
        #self.snap_number = snap_number
        self.outpath = outpath #"../plots/{}/".format(sim)
        self.analysis_type = analysis_type
        self.figure_name = figure_name

        self.simulation_directory = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100'.format(self.sim)

        # kwargs

        if 'remove_subs' in kwargs.keys():
            self.remove_subs = kwargs['remove_subs']
        else:
            self.remove_subs = False
        # Substructure

        if self.remove_subs == True :
            print("-> Removing particles from the massive subhalo using Emily Cunninghan's SCF IDs")
            subs_path = '/mnt/home/ecunningham/ceph/latte/{}_res7100/massive_stream/dm_inds.npy'.format(self.sim)
            subs_ids = np.load(subs_path)
            self.mask_sub=np.ones(len(dm_host_distance), dtype=bool) 
            self.mask_sub[subs_ids]=0
            self.subs_id =np.zeros_like(self.mask)
            self.subs_id[subs_ids]=1
        
    def analyze_snapshot(self, snap):  

        to_kpc_gyr = ((1*u.km/u.s).to(u.kpc/u.Gyr)).value
       
        # ** stellar particles ** 
        #stars_host_distance=part['star'].prop('host.distance.principal')
        #stars_host_velocity=part['star'].prop('host.velocity.principal')*to_kpc_gyr

        if  'star' in self.part_type:
            part = ga.io.Read.read_snapshots(self.part_type,'snapshot',
                snap, self.simulation_directory, assign_hosts=True,
                assign_orbits=True, assign_hosts_rotation=False)

            stars_host_distance=part['star'].prop('host.distance')
            stars_host_velocity=part['star'].prop('host.velocity')*to_kpc_gyr
            # Distance cuts
            stars_dist = np.sqrt(np.sum(stars_host_distance**2, axis=1))
            stars_dist_cut = np.where((stars_dist>self.rmin) & (stars_dist<self.rmax)) 
            pos_stars, vel_stars = rotate_halo(stars_host_distance[stars_dist_cut],
                                               stars_host_velocity[stars_dist_cut])
        

            if self.analysis_type == "orbital_poles":
                print("-> Computing orbital poles for star particles")
                stars_kinematics = nba.kinematics.Kinematics(pos_stars, vel_stars)
                stars_l_host, stars_b_host = stars_kinematics.pos_cartesian_to_galactic()
                stars_figname = self.outpath + "{}_".format(self.sim) + self.figure_name +  "_stars_{:03d}.png".format(snap)
                mollweide_projection(stars_l_host*180/np.pi, stars_b_host*180/np.pi, [0], [0], self.rmin, self.rmax,
                                     bmin=100, bmax=1000, nside=40, figname=stars_figname)
        
    
        # ** dark matter ** 
        if 'dark' in self.part_type:
            part = ga.io.Read.read_snapshots(self.part_type,'snapshot',
                snap, self.simulation_directory, assign_hosts=True,
                assign_orbits=True, assign_hosts_rotation=False,
                sort_dark_by_id=True)

            #dm_host_distance=part['dark'].prop('host.distance.principal')
            #dm_host_velocity=part['dark'].prop('host.velocity.principal')*to_kpc_gyr
        
            dm_host_distance=part['dark'].prop('host.distance')
            dm_host_velocity=part['dark'].prop('host.velocity')*to_kpc_gyr
	    
            #cartessian_projection(dm_host_distance, figure_name_dark+"with")
            # remove substructue. 
            if self.remove_subs == True :
                dm_host_distance = dm_host_distance[self.mask_sub]
                dm_host_velocity = dm_host_velocity[self.mask_sub]

                #m_subs_distance = np.mean(dm_host_distance[self.subs_id])
                #dm_subs_velocity = np.mean(dm_host_velocity[self.subs_id])

            #cartessian_projection(dm_host_distance[mask], figure_name_dark)
 
            dm_dist = np.sqrt(np.sum(dm_host_distance**2, axis=1))
            dm_dist_cut = np.where((dm_dist>self.rmin) & (dm_dist<self.rmax)) 
            pos_dm, vel_dm = rotate_halo(dm_host_distance[dm_dist_cut], dm_host_velocity[dm_dist_cut])

            if self.analysis_type == "orbital_poles":
                dm_m12b_kinematics_2 = nba.kinematics.Kinematics(pos_dm, vel_dm)
                dm_l_host, dm_b_host = dm_m12b_kinematics_2.pos_cartesian_to_galactic()
                dm_figname = self.outpath + "{}_".format(self.sim) + self.figure_name + "_{:03d}.png".format(snap)
                mollweide_projection(dm_l_host*180/np.pi, dm_b_host*180/np.pi,
                    [0],[0], self.rmin, self.rmax, bmin=10, bmax=1500,
                    nside=40, figname=dm_figname)
        return 0
    def callback(self, result):
        #ith open(self.outpath+"text.txt", 'a') as f:
        #self.write("{0}\n".format(result))
        print('here')

    def __call__(self, task):
        return self.analyze_snapshot(task)

def main(pool, sim, rmin, rmax, part_type, analysis_type, outpath, figure_name):
    worker = FIRE_analysis(sim, rmin, rmax, part_type, analysis_type, outpath, figure_name)
    tasks = np.arange(snap_init, snap_final, delta_snap)

    for r in pool.map(worker, tasks, callback=worker.callback):
        pass

    pool.close()


if __name__ == "__main__":
    #snap_final = 600
    #snap_init = 425
    #snap_final = 601
    delta_snap = 1
    #figure_name_dark = 'OP_dark_particles_shells_no_rot_no_substructure_'
    #figure_name_stars = 'OP_stars_particles_within_50kpc_'
    #outpath = "../plots/{}/".format(sim)
    #sim = 'm12b'

    #rmin = 50
    #rmax = 500
    part_type = ['dark']
    #snap_number = 385
    #figure_name = "OP_dark_particles_shells_no_rot_no_substructure_"
    outpath = "./" #./plots/{}/".format(sim)
    analysis_type = "orbital_poles"

    #m12b_OP = FIRE_analysis(sim, rmin, rmax, part_type, snap_number, analysis_type, outpath, figure_name)
    #m12b_OP.analyze_snapshot()
    

    from argparse import ArgumentParser
    parser = ArgumentParser(description="")

    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    
    #.add_argument("--mpi", dest="mpi", default=False,
    #                   action="store_true", help="Run with MPI.")
    
    parser.add_argument("--sim", dest="sim", default=False, 
                       help="simulation to use: m12b, m12c, m12i..")
    parser.add_argument("--figname", dest="figure_name", default=False, 
                       type=str, help="Simulation name")
    parser.add_argument("--rmin", dest="rmin", default=0, 
                       type=int, help="minimum radii")
    parser.add_argument("--rmax", dest="rmax", default=500, 
                       type=int, help="maximum radii")
    parser.add_argument("--i", dest="snap_init", default=0, 
                       type=int, help="final snapshot")
    parser.add_argument("--f", dest="snap_final", default=1, 
                       type=int, help="final snapshot")
    args = parser.parse_args()

   
    pool = schwimmbad.choose_pool(mpi=False, processes=args.n_cores)


    snap_init = args.snap_init
    snap_final = args.snap_final

    compute_poles = main(pool, args.sim, args.rmin, args.rmax,
                         part_type, analysis_type, outpath, args.figure_name,
                         remove_subs=True)
