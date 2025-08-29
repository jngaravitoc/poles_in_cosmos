import numpy as np
import matplotlib.pyplot as plt
import h5py
#import nba
import sys

sys.path.append("/mnt/home/ecunningham/python")
sys.path.append("../scripts/src/")
plt.style.use('~/matplotlib.mplstyle')

#import pynbody
import pynbody_routines as pr
import gizmo_analysis as ga

from matplotlib import colors
import halo_analysis as halo
import io_gizmo_pynbody as fa

from scipy.linalg import norm


# Get all the subhalos


# Get all the subhalos
def get_halo_satellite(sim, mass_rank, init_snap=300, final_snap=600):
    sim_directory = "/mnt/ceph/users/firesims/fire2/metaldiff/{}_res7100/".format(sim);
    m12_subhalos = halo.io.IO.read_catalogs('index', init_snap, sim_directory);

    halt = halo.io.IO.read_tree(simulation_directory=sim_directory);
    hsub = pr.pynbody_subhalos(m12_subhalos)
    sat_id = np.argsort(hsub.dark['mass'])[mass_rank]
    sat_tree_id = m12_subhalos['tree.index'][sat_id]
    satellite = fa.return_tracked_pos(halt, sat_tree_id, pynbody_halo=False, 
                                      init_snap=init_snap, final_snap=final_snap)
    return satellite


def write_group(cat_name, sat_id, subhalos):
    hf = h5py.File(cat_name, 'r+')
    snap_group = hf.create_group('{}'.format(sat_id))
    snap_group['mass'] = subhalos['mass']
    snap_group['pos']= subhalos['position']
    snap_group['vel'] = subhalos['velocity']
    snap_group['tree_index'] = subhalos['treeind']
    snap_group['snapshot'] = subhalos['snapshot']
    snap_group['catalogue_index'] = subhalos['catind']
    return 0

if __name__ == "__main__":
    sim = 'm12w'
    cat_name = '{}_satellites_catalog.h5'.format(sim)
    hf = h5py.File(cat_name, 'w')

    m12m_sat_list = []
    k=11
    for i in range(1, 300):
        m12m_sat1 = get_halo_satellite(sim, -i, init_snap=300, final_snap=600);
        if len(m12m_sat1['position']) > 0:
            m12m_sat_dist_z0 = norm(m12m_sat1['position'], axis=1)[-1]
            if m12m_sat_dist_z0 < 300:
                m12m_sat_list.append(i)
                write_group(cat_name, k, m12m_sat1)
                k+=1
    print('Done making satellite_catlog with a total of {} satellites'.format(len(m12m_sat_list)))
