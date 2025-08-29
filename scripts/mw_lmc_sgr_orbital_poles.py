import numpy as np
from nba.ios import load_halo, load_snapshot, get_com
from nba.visuals import Visuals
from nba.kinematics import Kinematics
import nba.com as com

if __name__ == "__main__":
    snapshot="/mnt/home/jhunt/ceph/data/Laporte/trial_mpa/snaps_lmc_sgr_h2/"
    snapname="snap"
    nhost=40000000
    nsat=1500000
    

    host_halo = load_halo(snapshot+snapname+"_001",  N_halo_part=[nhost, nsat], q=['pos', 'vel', 'mass'],  com_frame=0, galaxy=0,
        snapformat=1, com_method='shrinking')


    pos_com_host = host_halo['pos']
    vel_com_host = host_halo['vel']
    d_host = np.sqrt(np.sum(pos_com_host**2, axis=1))
    rcut = np.where( (d_host < 500) & (d_host > 50))

    halo_kinematics = Kinematics(pos_com_host[rcut], vel_com_host[rcut])
    op_l, op_b = halo_kinematics.orbpole()

    halo_vis = Visuals()
    twd_map, hpx_map_i = halo_vis.compute_mollweide(op_l, op_b, nside=48)
    
    for k in range(0, 101, 20):
        """
        pos_com_host, vel_com_host = load_halo(snapshot+snapname+"_{:03d}".format(k), N_halo_part=[nhost, nsat], q=['mass'],  com_frame=0, galaxy=0, snapformat=3, com_method='diskpot')
        d_host = np.sqrt(np.sum(pos_com_host**2, axis=1))
        rcut = np.where( (d_host < 500) & (d_host > 50))

        halo_kinematics = Kinematics(pos_com_host[rcut], vel_com_host[rcut])
        op_l, op_b = halo_kinematics.orbpole()
        halo_vis = Visuals()
        twd_map, hpx_map = halo_vis.compute_mollweide(op_l, op_b, nside=48)
        #twd_map = halo_vis.compute_mollweide(op_l, op_b, nside=68)
        figname = "OP_enhancement_MWLMC5_100M_b0_OM3_{:03d}.png".format(k)
        halo_vis.plot_mollweide_galactic(hpx_map/hpx_map -1, rotation=(180, 0, 0), bmin=-0.3, bmax=0.3, figname=figname, fig_title="")
        """


        """ 
        pos_disk = load_snapshot(snapshot+snapname+"_{:03d}".format(k), 1, 'pos', 'disk')
        vel_disk = load_snapshot(snapshot+snapname+"_{:03d}".format(k), 1, 'vel', 'disk')
        pos_com, vel_com = get_com(pos_disk, vel_disk, np.ones(len(pos_disk)), 'diskpot', snapshot+snapname+"_{:03d}".format(k), 3)
        new_pos = com.re_center(pos_disk, pos_com)
        new_vel = com.re_center(vel_disk, vel_com)
            
        disk_kinematics = Kinematics(new_pos, new_vel)

        op_l_d, op_b_d = disk_kinematics.orbpole()
        disk_vis = Visuals()
        twd_map, hpx_map = disk_vis.compute_mollweide(op_l_d, op_b_d, nside=48)
        #twd_map = halo_vis.compute_mollweide(op_l, op_b, nside=68)
        figname = "OP_center_disk_MWLMC5_100M_b0_OM3_{:03d}.png".format(k)
        disk_vis.plot_mollweide_galactic(hpx_map, rotation=(180, 0, 0),  bmin=0.0, bmax=100, figname=figname, fig_title="MWLMC_disk_OP")
        """
