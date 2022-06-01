import numpy as np
import utilities as ut
from satellite_analysis import halo_reader as hr
from satellite_analysis import satellite_io as sio


# m12 info
m12_names = ['m12f', 'm12b', 'm12c', 'm12w']

m12_dirs_peloton = ['/home/jsamuel/scratch/m12f/m12f_res7100',
                    '/home/jsamuel/scratch/m12b/m12b_res7100',
                    '/home/jsamuel/scratch/m12c/m12c_res7100',
                    '/home/jsamuel/scratch/m12w/m12w_res7100']


def track_main_host_properties(hal, hal_mask, host_str='host.'):
    host_ind = np.array([hal[host_str+'index'][hal_mask][0]])
    host_track_ = sio.halo_track(hal, host_ind, snapshot_list=np.arange(0,600,1))
    host_track = np.array([track[0] for track in host_track_], dtype='int32')

    host_history = {
        'snapshot':np.arange(1,601,1),
        'host.index':host_track,
        'host.mass':np.array(hal['mass'][host_track]),
        'host.radius':np.array(hal['radius'][host_track]),
        'host.scale.radius':np.array(hal['scale.radius'][host_track]),
        'host.position':hal['position'][host_track],
        'host.velocity':hal['velocity'][host_track]
    }

    return host_history

# satellite selection is arbitrary since I'm only interested in the host here
m12_01 = hr.SatelliteTree(directory_list=m12_dirs_peloton,
                        host_name_list=m12_names,
                        mask_names=['most.star.mass'],
                        assign_species=True,
                        snapshot_indices=[600],
                        number_sats=14)

m12_host_history = sio.loop_hal(m12_01, 'most.star.mass', track_main_host_properties)

# save host halo history separately
for host in m12_names:
    host_history = {}
    for history_key in m12_host_history[host][0].keys():
        if history_key.rsplit('.')[0] == 'host':
            host_history[history_key] = m12_host_history[host][0][history_key]
    host_history['snapshot'] = np.arange(0,601,1,dtype=np.int32)
    host_history['time'] = np.array(
        sio.convert_snapshot_to_time(
            np.arange(0,601,1), 
            time_kind='time', 
            snap_file_path='/home/jsamuel/tables/'),
        dtype=float)
    host_history['redshift'] = np.array(
        sio.convert_snapshot_to_time(
            np.arange(0,601,1), 
            time_kind='redshift', 
            snap_file_path='/home/jsamuel/tables/'),
        dtype=float)
    ut.io.file_hdf5('{}_host_halo_history'.format(host), dict_or_array_to_write=host_history)
