import numpy as np
from nba.ios import load_halo, load_snapshot, get_com
from nba.visuals import Visuals
from nba.kinematics import Kinematics
import nba.com as com

if __name__ == "__main__":
    #snapshot = "/mnt/home/nico/ceph/gadget_runs/MWLMC/MWLMC5_b0/out/"
    #snapname = 'MWLMC5_100M_b0_vir_OM3_G4'
    snapshot="/mnt/home/nico/ceph/MWLMC_sims/MWLMC/MWLMC5_b0/"
    snapname="MWLMC5_100M_beta0_vir_OM5_G4"
    nhost=100000000-1
    nsat=15000000
    k=000

    pos_disk = load_snapshot(snapshot+snapname+"_{:03d}".format(k), 3, 'pos', 'disk')
    acc_disk = load_snapshot(snapshot+snapname+"_{:03d}".format(k), 3, 'acc', 'dm')
    vel_disk = load_snapshot(snapshot+snapname+"_{:03d}".format(k), 3, 'vel', 'disk')
    pos_com, vel_com = get_com(pos_disk, vel_disk, np.ones(len(pos_disk)), method='shrinking')
    new_pos = com.re_center(pos_disk, pos_com)

    print(acc_disk[0])
    d_host = np.sqrt(np.sum(new_pos**2, axis=1))
    a_host = np.sqrt(np.sum(acc_disk**2, axis=1))
    rcut = np.where( (d_host < 5) )


    fig = plt.figure()
    plt.plot(d_host[rcut], a_host[rcut] )
    plt.savefig('a_profile.png')
    plt.close()
    print('Done')
