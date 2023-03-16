"""
Routines to use FIRE simulations
author: jngaravitoc, with contributions from Emily Cunningham, Arpit Arora
github: jngaravitoc

02/2022 - 



TODO:
-----
- Select stars of the disk using Emily's routine [Done] ->  cant his only be
  used in snap 600?

"""


import numpy as np
from numpy import Inf

import sys
import schwimmbad
from scipy.optimize import curve_fit
from astropy import units as u
#import gala.potential as gp
#from gala.units import galactic, solarsystem, dimensionless
#from gala.potential import scf

from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys
sys.path.append("/mnt/home/ecunningham/python")
sys.path.append("/mnt/home/nico/codes/nbody_analysis/src")

import gizmo_analysis as ga
import halo_analysis as ra
import utilities as ut
import nba
from scipy.linalg import norm

import pynbody
from pynbody import transformation
import pynbody_routines as pr

## Tracking subhalos using their index at 300th snap (using only merger tree)
def return_tracked_pos(halo_tree, tr_ind_at300, pynbody_halo=False):
    # Adapted from Arpit's function
    h_index = tr_ind_at300
    tree_ind = []
    for _ in range(300,601):
        tree_ind.append(h_index)
        h_index = halo_tree['descendant.index'][h_index]
    tree_ind = np.array(tree_ind)
    position = halo_tree['host.distance'][tree_ind]
    nsnaps = halo_tree['snapshot'][tree_ind]
    mass = halo_tree['mass'][tree_ind]
    velocity = halo_tree['host.velocity'][tree_ind]
    #vel_rad = halt['host.velocity.rad'][tree_ind]
    #vel_tan = halt['host.velocity.tan'][tree_ind]
    if pynbody_halo == True:
        sat = {'position': position,
               'mass': mass,
               'velocity': velocity,
            }
        return pr.pynbody_satellite(sat)

    elif pynbody_halo == False:

      return {'position': position,
              'snaps' : nsnaps,
              'mass' : mass,
              'velocity' : velocity,
             }


class FIRE_analysis(object):
    def __init__(self, sim, rmin, rmax, part_type, analysis_type, outpath, figure_name, remove_subs, **kwargs):
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
        self.outpath = outpath #"../plots/{}/".format(sim)
        self.analysis_type = analysis_type
        self.figure_name = figure_name

        self.simulation_directory = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100'.format(self.sim)
        times = '/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt'.format(self.sim)
        self.times = np.loadtxt(times, usecols=3)
        self.remove_subs = remove_subs
        print(self.remove_subs)
         
        self.bmin = kwargs['bmin'] 
        self.bmax = kwargs['bmax']

        if self.remove_subs == 1 :
            print("-> Removing particles from the massive subhalo using Emily Cunninghan's IDs")
            subs_path = '/mnt/home/ecunningham/ceph/latte/{}_res7100/massive_stream/dm_inds.npy'.format(self.sim)
            self.subs_ids = np.load(subs_path)
            #self.mask_sub=np.ones(len(dm_host_distance), dtype=bool) 
            #self.mask_sub[subs_ids]=0
            #self.subs_id =np.zeros_like(self.mask)
            #self.subs_id[subs_ids]=1
        
    def analyze_snapshot(self, snap):  

        to_kpc_gyr = ((1*u.km/u.s).to(u.kpc/u.Gyr)).value
       
        # ** stellar particles ** 
        #stars_host_distance=part['star'].prop('host.distance.principal')
        #stars_host_velocity=part['star'].prop('host.velocity.principal')*to_kpc_gyr

        if 'disk' in self.part_type:
            part = ga.io.Read.read_snapshots(['dark', 'star'],'snapshot',
                snap, self.simulation_directory, assign_hosts=True,
                assign_orbits=True)

            stars_host_distance=part['star'].prop('host.distance')
            stars_host_velocity=part['star'].prop('host.velocity')*to_kpc_gyr
            dm_host_distance=part['dark'].prop('host.distance')
            dm_host_velocity=part['dark'].prop('host.velocity')*to_kpc_gyr

            dm_host_mass=part['dark'].prop('mass')
            stars_host_mass=part['star'].prop('mass')

            #--- Disk stars

            
            #form_dist=part['star'].prop('form.host.distance.total')
            #form_time=part['star'].prop('form.time')
            disk_particles = np.where(stars_host_distance<30)[0]
            #print(len(disk_particles))
            pos_disk = stars_host_distance[disk_particles]
            vel_disk = stars_host_velocity[disk_particles]
            disk_mass = stars_host_mass[disk_particles]

            #--- pynbody 
            
            f = pynbody_halo(stars_host_distance, stars_host_velocity,
                stars_host_mass, dm_host_distance, dm_host_velocity,
                dm_host_mass, pos_disk, vel_disk, disk_mass)
            path="/mnt/home/nico/projects/poles_in_cosmos/scripts/"
            #pynbody.plot.hist2d(f.stars[:,0], f.stars[:,1], make_plot=True,
            #    cmap='Blues', filename=path+'pynbody_hist_plot.png')
            print("here pynbody") 
            center = pynbody.analysis.halo.center(f, mode='ssc')
            #print("center", transformation.inverse_translate(center, [39762.54, 41743.92, 39791.35]))
            #enter2 = pynbody.analysis.halo.center(f, mode='ssc', retcen=True)
            #rint("center", center2)  
            
            center3 = pynbody.analysis.angmom.sideon(f, cen=(0,0,0))
            ## todo -> load halo centers and make the first element  the halo
            ## center and then revert and print the zeroth element. 
            sim_center = np.loadtxt("m12b_center.txt")
            f['pos'][0] = sim_center[snap]
            print("here", sim_center[snap])
            center3.revert()
            halo_center_xy_plane = f['pos'][0]
            #print(f['pos'][0])
            
            """
            cartessian_projection(f.gas['pos'], "m12b_disk_pynbody_centered{:03d}.png".format(snap))
            disk_kinematics = nba.kinematics.Kinematics(f.gas['pos'], f.gas['vel'])
            disk_l_host, disk_b_host = disk_kinematics.pos_cartesian_to_galactic()
            OP_disk_l_host, OP_disk_b_host = disk_kinematics.orbpole()
            mollweide_projection(disk_l_host*180/np.pi, disk_b_host*180/np.pi, [0], [0], 
                                 title="disk particles", bmin=30, bmax=1000,
                                 nside=40, figname="m12b_mollweide_disk_pynbody_centered_{:03d}.png".format(snap))

           
            mollweide_projection(OP_disk_l_host, OP_disk_b_host, [0], [0], 
                                 title="disk orbital poles", bmin=30, bmax=1000,
                                 nside=40, figname="m12b_OP_disk_pynbody_centered_{:03d}.png".format(snap))

            disk_kinematics_off = nba.kinematics.Kinematics(pos_disk, vel_disk)
            disk_l_host_off, disk_b_host_off = disk_kinematics_off.pos_cartesian_to_galactic()
            OP_disk_l_host_off, OP_disk_b_host_off = disk_kinematics_off.orbpole()
            mollweide_projection(disk_l_host_off*180/np.pi, disk_b_host_off*180/np.pi, [0], [0], 
                                 title="disk particles", bmin=30, bmax=1000,
                                 nside=40,  figname="m12b_mollweide_disk_pynbody_off_centered_{:03d}.png".format(snap))

            mollweide_projection(OP_disk_l_host_off, OP_disk_b_host_off, [0], [0], 
                                 title="disk orbital poles", bmin=30, bmax=1000,
                                 nside=40,  figname="m12b_OP_disk_pynbody_off_centered_{:03d}.png".format(snap))
@           """
        if 'star' in self.part_type:
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
                stars_figname = self.outpath + "{}_".format(self.sim) + self.figure_name +  "_{:03d}.png".format(snap)
                title_stars = "{} star; {}-{} kpc; t={:.2f}  Gyr".format(self.sim, self.rmin,  self.rmax, self.times[snap] )
                mollweide_projection(stars_l_host*180/np.pi, stars_b_host*180/np.pi, [0], [0], 
                                     title=title_stars, bmin=self.bmin, bmax=self.bmax, nside=40, figname=stars_figname)
        
       
        # ** dark matter ** 
        if 'dark' in self.part_type:
            part = ga.io.Read.read_snapshots(['dark', 'star'],'snapshot',
                snap, self.simulation_directory, assign_hosts=True,
                assign_orbits=True, assign_hosts_rotation=False,
                sort_dark_by_id=True)

            #dm_host_distance=part['dark'].prop('host.distance.principal')
            #dm_host_velocity=part['dark'].prop('host.velocity.principal')*to_kpc_gyr
        
            stars_host_distance=part['star'].prop('host.distance')
            stars_host_velocity=part['star'].prop('host.velocity')*to_kpc_gyr
            dm_host_distance=part['dark'].prop('host.distance')
            dm_host_velocity=part['dark'].prop('host.velocity')*to_kpc_gyr

            dm_host_mass=part['dark'].prop('mass')
            stars_host_mass=part['star'].prop('mass')
            
            disk_particles = np.where(stars_host_distance<30)[0]
            #print(len(disk_particles))
            pos_disk = stars_host_distance[disk_particles]
            vel_disk = stars_host_velocity[disk_particles]
            disk_mass = stars_host_mass[disk_particles]
            
            f = pynbody_halo(stars_host_distance, stars_host_velocity,
                stars_host_mass, dm_host_distance, dm_host_velocity,
                dm_host_mass, pos_disk, vel_disk, disk_mass)

            path="/mnt/home/nico/projects/poles_in_cosmos/scripts/"
            #pynbody.plot.hist2d(f.stars[:,0], f.stars[:,1], make_plot=True,
            #    cmap='Blues', filename=path+'pynbody_hist_plot.png')
            center = pynbody.analysis.halo.center(f, mode='ssc')
            pynbody.analysis.angmom.faceon(f, cen=(0,0,0))
	    
            #cartessian_projection(dm_host_distance, figure_name_dark+"with")
            # remove substructue. 
            if self.remove_subs == 1 :
                mask_sub=np.ones(len(dm_host_distance), dtype=bool) 
                mask_sub[self.subs_ids]=0
                subs_id = np.zeros_like(mask_sub)
                subs_id[self.subs_ids]=1

                dm_host_distance = dm_host_distance[mask_sub]
                dm_host_velocity = dm_host_velocity[mask_sub]

                #m_subs_distance = np.mean(dm_host_distance[self.subs_id])
                #dm_subs_velocity = np.mean(dm_host_velocity[self.subs_id])

            #cartessian_projection(dm_host_distance[mask], figure_name_dark)
 
            dm_dist = np.sqrt(np.sum(dm_host_distance**2, axis=1))
            dm_dist_cut = np.where((dm_dist>self.rmin) & (dm_dist<self.rmax)) 
            #os_dm, vel_dm = rotate_halo(dm_host_distance[dm_dist_cut], dm_host_velocity[dm_dist_cut])

            if self.analysis_type == "orbital_poles":
                dm_m12b_kinematics_2 = nba.kinematics.Kinematics(dm_host_distance[dm_dist_cut], dm_host_velocity[dm_dist_cut])
                dm_l_host, dm_b_host = dm_m12b_kinematics_2.pos_cartesian_to_galactic()
                OP_dm_l_host, OP_dm_b_host = dm_m12b_kinematics_2.orbpole()
                dm_figname = self.outpath + "{}_".format(self.sim) + self.figure_name + "_{:03d}.png".format(snap)
                dm_figname2 = self.outpath + "OP_{}_".format(self.sim) + self.figure_name + "_{:03d}.png".format(snap)
                title_dark = "{} dark; {}-{} kpc; t={:.2f} Gyr".format(self.sim, self.rmin, self.rmax, self.times[snap])
                mollweide_projection(dm_l_host*180/np.pi, dm_b_host*180/np.pi,
                    [0],[0], title=title_dark, bmin=self.bmin,  bmax=self.bmax,
                    nside=40, figname=dm_figname)

                mollweide_projection(OP_dm_l_host, OP_dm_b_host,
                    [0],[0], title=title_dark, bmin=self.bmin,  bmax=self.bmax,
                    nside=40, figname=dm_figname2)

            if self.analysis_type == "shape":
                shape = pynbody.analysis.halo.halo_shape(pos_dm, N=100, rin=20, rout=500, bins='equal')
                #dm_figname = self.outpath + "{}_".format(self.sim) + self.figure_name + "_{:03d}.png".format(snap)
                title_dark = "{} dark {}-{} kpc".format(self.sim, self.rmin,  self.rmax)
        
        return halo_center_xy_plane
        
    def callback(self, result):
        with open(self.outpath+"text.txt", 'a') as f:
            f.write("{0}\n".format(result))
        #print('here')

    def __call__(self, task):
        return self.analyze_snapshot(task)

def main(pool, sim, rmin, rmax, part_type, analysis_type, outpath, figure_name,  remove_subs, **kwargs): 
    bmin = kwargs['bmin']
    bmax = kwargs['bmax']

    worker = FIRE_analysis(sim, rmin, rmax, part_type, analysis_type, outpath, figure_name, remove_subs, bmin=bmin, bmax=bmax)
    print(worker)
    tasks = np.arange(snap_init, snap_final, delta_snap)

    for r in pool.map(worker, tasks, callback=worker.callback):
        pass

    pool.close()

def write_params(params):
    """
    """
    print("-> Writing parameter file in ")
    
    filename = params.figure_name
    dt = datetime.now()    
    outpath = params.outpath
    f = open(outpath+"parameter_file_" + filename + ".txt", "w")
    f.write(str(dt) + "+ \n")
    f.write(str(params))
    f.close()




if __name__ == "__main__":
    #snap_final = 600
    #snap_init = 425
    #snap_final = 601
    #figure_name_dark = 'OP_dark_particles_shells_no_rot_no_substructure_'
    #figure_name_stars = 'OP_stars_particles_within_50kpc_'
    #outpath = "../plots/{}/".format(sim)
    #sim = 'm12b'
    #rmin = 50
    #rmax = 500
    #snap_number = 385
    #figure_name = "OP_dark_particles_shells_no_rot_no_substructure_"
    #outpath = "./" #./plots/{}/".format(sim)

    delta_snap = 1
    analysis_type = "orbital_poles"

    #m12b_OP = FIRE_analysis(sim, rmin, rmax, part_type, snap_number, analysis_type, outpath, figure_name)
    #m12b_OP.analyze_snapshot()
    

    from argparse import ArgumentParser
    parser = ArgumentParser(description="")

    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    
    #roup.add_argument("--mpi", dest="mpi", default=False,
    #                   action="store_true", help="Run with MPI.")
    
    parser.add_argument("--sim", dest="sim", default=False, 
                       help="simulation to use: m12b, m12c, m12i..")
    parser.add_argument("--figname", dest="figure_name", default=False, 
                       type=str, help="Simulation name")
    parser.add_argument("--outpath", dest="outpath", default="./", 
                       type=str, help="output path")
    parser.add_argument("--rmin", dest="rmin", default=0, 
                       type=int, help="minimum radii")
    parser.add_argument("--rmax", dest="rmax", default=500, 
                       type=int, help="maximum radii")
    parser.add_argument("--bmin", dest="bmin", default=10, 
                       type=int, help="minimum b for histogram")
    parser.add_argument("--bmax", dest="bmax", default=1500, 
                       type=int, help="maximum b for histrogram")
    parser.add_argument("--i", dest="snap_init", default=0, 
                       type=int, help="final snapshot")
    parser.add_argument("--f", dest="snap_final", default=1, 
                       type=int, help="final snapshot")
    parser.add_argument("--remove_subs", dest="remove_subs", 
                       type=int, default=0,  help="removing substructure")
    parser.add_argument("--partype", dest="part_type", 
                       type=int, help="particle type 1) DM, 2) stars, 3) both")
    args = parser.parse_args()

   
    pool = schwimmbad.choose_pool(mpi=False, processes=args.n_cores)


    snap_init = args.snap_init
    snap_final = args.snap_final



    if args.part_type == 1:
       part_type=['dark']
       #args.bmin=5
       #args.bmax=7
       #args.bmax=9
    elif args.part_type == 2:
       #part_type=['star']
       part_type=['disk']
       #args.bmin=0
       #args.bmax=5

    elif args.part_type == 3:
       part_type=['dark', 'star']
    
    write_params(args) 

    # Other parameters 
    # bmin, bmax = 100, 1500 -> dark
    # bmin, bmax = 10, 100 -> stars
    

    delta_snap = 1
    analysis_type = "orbital_poles"

    halo_center_xy_plane = np.zeros(601)
    compute_poles = main(pool, args.sim, args.rmin, args.rmax, part_type,
                         analysis_type, args.outpath, args.figure_name,
                         args.remove_subs, bmin=args.bmin, bmax=args.bmax)

