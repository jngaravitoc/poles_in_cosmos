import numpy as np
import sys

sys.path.append("/mnt/home/ecunningham/python")

import io_sims as ios



sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/m12b_res7100/"
ios.read_snapshot(simname="FIRE", directory=sim_directory, snapnum=385, 
                  partypes={
                      'species' : ['dark'],
                      'snapshot_values_kind' : 'snapshots',
                      'snapshot_values' : 385,
                      'simulation_directory' : sim_directory,
                      'snapshot_directory' : '',
                      'simulation_name' : 'm12b',
                      'properties' : 'all',
                      'particle_subsample_factor' : None,
                      'sort_dark_by_id' : False,
                      'host_number' : 1,
                      'assign_hosts' : True,
                      'assign_orbits' : False,
                      'assign_formation_coordinates' : False,
                      
                  }
                  
                 )





