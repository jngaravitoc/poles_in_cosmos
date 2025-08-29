"""

Author: Nico Garavito-Camargo
Github: jngaravitoc

"""
#!/usr/bin/env python
# coding: utf-8


import numpy as np
import sys
from io_gizmo_pynbody import FIRE

sys.path.append("/mnt/home/ecunningham/python")


if __name__ == "__main__":
    
    # Set parameters 

    snap_init = 386
    snap_final = 387
    sim='m12b'
    
    
    # Data paths 
   

    snap_times = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/snapshot_times.txt".format(sim)
    times = np.loadtxt(snap_times, usecols=3)
    

    #data = FIRE(sim, remove_satellite=True, rm_stellar_sat=True)
    data_sat = FIRE(sim, only_sat=True)
    


    for k in range(snap_init, snap_final, 1):
        print(times[k])
        #hfaceon, hedge_on_ = data.rotated_halo(snap=k)
        hfaceon, hedge_on_ = data_sat.rotated_halo(snap=k)
        #pos_dm = hfaceon.dark['pos']
        #vel_dm = hfaceon.dark['vel']
        pos_dm_sat = hfaceon.dark['pos']
        vel_dm_sat = hfaceon.dark['vel']
        #pos_st = hfaceon.stars['pos']
        #vel_st = hfaceon.stars['vel']
        #data_dm = np.array([pos_dm[:,0], pos_dm[:,1], pos_dm[:,2], vel_dm[:,0], vel_dm[:,1], vel_dm[:,2]]) 
        data_dm_sat = np.array([pos_dm_sat[:,0], pos_dm_sat[:,1], pos_dm_sat[:,2], vel_dm_sat[:,0], vel_dm_sat[:,1], vel_dm_sat[:,2]]) 
        #data_st = np.array([pos_st[:,0], pos_st[:,1], pos_st[:,2], vel_st[:,0], vel_st[:,1], vel_st[:,2]])
        #print(np.shape(data_dm), np.shape(data_st))
        np.save('dm_m12b_satellite_particles_snap_{:03d}'.format(k), data_dm_sat.T)
        #np.save('dm_m12b_halo_particles_snap_{:03d}'.format(k), data_dm.T)
        #np.save('stars_m12b_halo_particles_snap_{:03d}'.format(k), data_st.T)
