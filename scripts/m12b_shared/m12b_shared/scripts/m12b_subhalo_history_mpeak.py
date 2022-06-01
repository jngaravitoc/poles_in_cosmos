import numpy as np
from collections import defaultdict
from numba import jit
import utilities as ut
from satellite_analysis import halo_reader as hr
from satellite_analysis import satellite_io as sio


def halo_track(tree, initial_tree_ids, redshift_list=None, snapshot_list=None):
    """
    Track halos in a merger tree at specific times given by redshift or
    given by snapshot.
    """
    if redshift_list:
        snapshot_ids = tree.Snapshot.get_snapshot_indices('redshift', redshift_list)
    else:
        snapshot_ids = snapshot_list

    hal_tracker = np.full((601, len(initial_tree_ids)), -1, dtype='int32')
    i = np.max(snapshot_ids)
    hal_tracker[i] = initial_tree_ids

    tracking_ids = initial_tree_ids
    while i >= np.min(snapshot_ids):
        i -= 1
        progenitor_indices = tree['progenitor.main.index'][tracking_ids]
        if i in snapshot_ids:
            tracking_ids = optim_track(progenitor_indices)
            hal_tracker[i] = tracking_ids

    return hal_tracker

@jit(nopython=True)
def optim_track(progenitor_inds):
    negative_ids = np.where(progenitor_inds < 0)[0]
    progenitor_inds[negative_ids] = -1#this may not work for all trees?

    return progenitor_inds


def infall_hist(hal, hal_mask, host_str='host.', min_snapshot=0, time_file_path='.'):
    # hal_mask should be for a selection of satellites at a single snapshot
    initial_snapshot = hal['snapshot'][hal_mask][0]
    all_snapshots = np.arange(min_snapshot, initial_snapshot+1, 1)
    initial_ids = np.where(hal_mask)[0]# tree indices

    # track the selection of satellites from z=0 all the way back to min_snapshot
    # return their tree indices at all snapshots
    tracked_ids = halo_track(hal, initial_ids, snapshot_list=all_snapshots)

    # also track the main MW host halo and return its tree indices at all
    # snapshots
    host_ind = np.array([hal[host_str+'index'][hal_mask][0]])
    host_track_ = sio.halo_track(hal, host_ind, snapshot_list=all_snapshots)
    host_track = np.array([track[0] for track in host_track_], dtype='int32')

    # get the central index for each tracked satellite at all snapshots
    central_mask = np.array(hal['central.index'][tracked_ids]) >= 0
    central_ind = np.where(
        central_mask, 
        np.array(hal['central.index'][tracked_ids]), 
        np.full(np.array(hal['central.index'][tracked_ids]).shape, -1))
    
    # set up nan arrays and 3D mask for things that don't have a central
    nan_central = np.full(np.array(tracked_ids).shape, np.nan)
    central_mask_3d = np.reshape(
        np.dstack([central_mask, central_mask, central_mask]), np.array(tracked_ids).shape + (3,))
    nan_central_3d = np.full(np.array(tracked_ids).shape + (3,), np.nan)


    infall_history = {
        'snapshot':np.arange(0,601,1),
        'time':np.array(
            sio.convert_snapshot_to_time(np.arange(0,601,1), 
                time_kind='time', snap_file_path=time_file_path), dtype=np.float),
        'redshift':np.array(
            sio.convert_snapshot_to_time(np.arange(0,601,1), 
                time_kind='redshift', snap_file_path=time_file_path), dtype=np.float),
        'tree.index':tracked_ids,
        'mass.200m':np.array(hal['mass'][tracked_ids]),
        'position':hal['position'][tracked_ids],
        'velocity':hal['velocity'][tracked_ids],
        'radius.200m':hal['radius'][tracked_ids],
        'nfw.scale.radius':hal['scale.radius'][tracked_ids],
        'vel.circ.max':hal['vel.circ.max'][tracked_ids],
        'star.mass':np.array(hal['star.mass'][tracked_ids]),
        'central.index':np.array(hal['central.index'][tracked_ids]),
        'central.mass.200m':np.where(central_mask, hal['mass'][central_ind], nan_central),
        'central.position':np.where(central_mask_3d, hal['position'][central_ind], nan_central_3d),
        'central.velocity':np.where(central_mask_3d, hal['velocity'][central_ind], nan_central_3d),
        'central.radius.200m':np.where(central_mask, hal['radius'][central_ind], nan_central),
        'central.nfw.scale.radius':np.where(central_mask, hal['scale.radius'][central_ind], nan_central),
        'central.vel.circ.max':np.where(central_mask, hal['vel.circ.max'][central_ind], nan_central),
        'central.star.mass':np.where(central_mask, hal['star.mass'][central_ind], nan_central)
    }

    return infall_history


def get_unique_subhalos(subhalo_history):
    unique_subhalos_dict = defaultdict(np.array)
    snapshot_array = np.arange(0, 601, 1)[::-1]
    subhalo_history = subhalo_history[::-1]

    for i,current_subhalo_track in enumerate(subhalo_history):
        if i == 0:
            for key in current_subhalo_track.keys():
                unique_subhalos_dict[key] = current_subhalo_track[key]
        else:
            current_tree_indices = current_subhalo_track['tree.index']

            j_mask = []
            for subhalo in current_tree_indices[snapshot_array[i]]:
                if (subhalo != -1) & (subhalo in unique_subhalos_dict['tree.index'].flatten()):
                    j_mask.append(False)
                else:
                    j_mask.append(True)

            if np.sum(j_mask) == 0:
                pass
            else:
                for key in current_subhalo_track.keys():
                    if key in ['snapshot', 'redshift', 'time']:
                        pass
                    else:
                        unique_subhalos_dict[key] = np.append(
                            unique_subhalos_dict[key], np.compress(
                                j_mask, current_subhalo_track[key], axis=1), axis=1)
                print(i, 'appended unique dictionary entries')

    lost_subhalos = unique_subhalos_dict['tree.index'] < 0
    lost_subhalos_3d = np.dstack([lost_subhalos, lost_subhalos, lost_subhalos])
    for key in unique_subhalos_dict.keys():
        if key in ['snapshot', 'redshift', 'time']:
            pass
        #elif key in ['position', 'velocity']:
        elif ('position' in key) | ('velocity' in key):
            unique_subhalos_dict[key] = np.where(
                lost_subhalos_3d, np.nan, unique_subhalos_dict[key])
        else:
            unique_subhalos_dict[key] = np.where(
                lost_subhalos, np.nan, unique_subhalos_dict[key])

    return unique_subhalos_dict


################################################################################

# Stampede2 locations
#m12b_dir = ['/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/m12b_res7100']
#time_tbl_path = '/home/jsamuel/tables/snapshot_times.txt'
# Peloton locations
m12b_dir = ['/home/jsamuel/scratch/m12b/m12b_res7100']
time_tbl_path = '/home/jsamuel/tables/'

### load data

m12_st = hr.SatelliteTree(directory_list=m12b_dir,
                            host_name_list=['m12b'],
                            mask_names=['mass.peak'],
                            assign_species=True,
                            mass_peak=1e8,
                            radius_limit=500,
                            redshift_limits=[0,1],
                            time_info_file_path=time_tbl_path+'snapshot_times.txt')

m12_subhalo_history = sio.loop_hal(m12_st, 'mass.peak', infall_hist, **{'time_file_path':time_tbl_path})

clean_m12_subhalo_history = get_unique_subhalos(m12_subhalo_history['m12b'])

ut.io.file_hdf5('m12b_subhalo_history_mpeak1e8_50-500kpc', dict_or_array_to_write=clean_m12_subhalo_history)

