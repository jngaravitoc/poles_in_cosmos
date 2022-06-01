import numpy as np
import utilities as ut
import gizmo_analysis as gizmo
from satellite_analysis import satellite_io as sio


m12_names = ['m12f', 'm12b', 'm12c', 'm12w']

m12_dirs_peloton = ['/home/jsamuel/scratch/m12f/m12f_res7100',
                    '/home/jsamuel/scratch/m12b/m12b_res7100',
                    '/home/jsamuel/scratch/m12c/m12c_res7100',
                    '/home/jsamuel/scratch/m12w/m12w_res7100']

for host, sim_dir in zip(m12_names, m12_dirs_peloton):
    part = gizmo.io.Read.read_snapshots(
        'star', 'snapshot', 600, simulation_directory=sim_dir,
        assign_hosts_rotation=True)

    disc_props_hdf5 = {}
    for key in ['position', 'velocity', 'rotation']:
        disc_props_hdf5[key] = part.hostz[key]

    disc_props_hdf5['snapshot'] = np.arange(0,601,1,dtype=np.int32)
    disc_props_hdf5['time'] = np.array(
        sio.convert_snapshot_to_time(
            np.arange(0,601,1), 
            time_kind='time', 
            snap_file_path='/home/jsamuel/tables/'),
        dtype=float)
    disc_props_hdf5['redshift'] = np.array(
        sio.convert_snapshot_to_time(
            np.arange(0,601,1), 
            time_kind='redshift', 
            snap_file_path='/home/jsamuel/tables/'),
        dtype=float)

    ut.io.file_hdf5('{}_host_disc_history'.format(host), dict_or_array_to_write=disc_props_hdf5)
