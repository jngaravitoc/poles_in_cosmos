import numpy as np
import sys
import nba
import pynbody


# Load sims snapshots

def load_halo(snapname, Npart_halo, galaxy):
    
    pos_dm = nba.ios.load_snapshot(snapname, snapformat=3, quantity='pos', ptype='dm')
    vel_dm = nba.ios.load_snapshot(snapname, snapformat=3, quantity='vel', ptype='dm')
    mass_dm = nba.ios.load_snapshot(snapname, snapformat=3, quantity='mass', ptype='dm')
    pids_dm = nba.ios.load_snapshot(snapname, snapformat=3, quantity='pid', ptype='dm')
    galaxy_ids = nba.ios.io_snaps.halo_ids(pids_dm, Npart_halo, galaxy)
    
    return pos_dm[galaxy_ids], vel_dm[galaxy_ids], mass_dm[galaxy_ids], pids_dm[galaxy_ids]

def load_disk(snapname):
    
    pos_disk = nba.ios.load_snapshot(snapname, snapformat=3, quantity='pos', ptype='bulge')
    vel_disk = nba.ios.load_snapshot(snapname, snapformat=3, quantity='vel', ptype='bulge')
    mass_disk = nba.ios.load_snapshot(snapname, snapformat=3, quantity='mass', ptype='bulge')
    pot_disk = nba.ios.load_snapshot(snapname, snapformat=3, quantity='pot', ptype='bulge')
    
    return pos_disk, vel_disk, mass_disk, pot_disk
    
    
def pynbody_satellite(pos, vel, mass):
    npart = len(mass)
    pynbody_sat = pynbody.new(dark=npart, order='dark')
    pynbody_sat.dark['pos'] = pos
    pynbody_sat.dark['vel'] = vel
    pynbody_sat.dark['mass'] = mass *1E10
    pynbody_sat.dark['pos'].units = 'kpc'
    pynbody_sat.dark['vel'].units = 'km s**-1'
    pynbody_sat.dark['mass'].units = 'Msol'
    return pynbody_sat
    
    
def pynbody_bulge(pos, vel, mass, pot):
    npart = len(mass)
    pynbody_bulge = pynbody.new(star=npart, order='star')
    pynbody_bulge.star['pos'] = pos
    pynbody_bulge.star['vel'] = vel
    pynbody_bulge.star['mass'] = mass
    pynbody_bulge.star['phi'] = pot

    pynbody_bulge.star['pos'].units = 'kpc'
    pynbody_bulge.star['vel'].units = 'km s**-1'
    pynbody_bulge.star['mass'].units = 'Msol'
    pynbody_bulge.star['phi'].units = 'km**2 s**-2'
    return pynbody_bulge


if __name__ == "__main__":
    sim = sys.argv[1]
    com_sat = np.zeros((80, 3))
    com_disk = np.zeros((80, 3))

    vel_sat = np.zeros((80, 3))
    vel_disk = np.zeros((80, 3))



    if sim == "LMC5":
        nsat = 15000000
    elif sim == "LMC6":
        nsat = 20830000
    elif sim =="LMC4":
        nsat = 8333333
    elif sim=="LMC3":
        nsat = 6666666


    j=0
    out_name_sat = '{}_orbit_satellite.txt'.format(sim)
    out_name_host = '{}_orbit_host.txt'.format(sim)

    for i in range(0, 400, 5):
        lmc5_path = '/mnt/home/nico/ceph/gadget_runs/MWLMC/MW{}_b0/out/MW{}_100M_b0_vir_OM3_G4_{:03d}'.format(sim, sim, i)
        sat = load_halo(lmc5_path, [100000000, nsat], 1)
        disk = load_disk(lmc5_path)
        pyn_sat = pynbody_satellite(sat[0], sat[1], sat[2])
        pyn_disk = pynbody_bulge(disk[0], disk[1], disk[2], disk[3])
        com_sat[j] = pynbody.analysis.halo.center(pyn_sat, mode='ssc',retcen=True)[:]
        com_disk[j] = pynbody.analysis.halo.center(pyn_disk, mode='pot',retcen=True)[:]
        pyn_disk.star['pos'] -= com_disk[j]
        pyn_sat.dark['pos'] -= com_sat[j]
        vel_disk[j] = pynbody.analysis.halo.vel_center(pyn_disk, retcen=True)
        vel_sat[j] = pynbody.analysis.halo.vel_center(pyn_sat, retcen=True)
        j+=1

    np.savetxt(out_name_sat, np.array([com_sat[:,0], com_sat[:,1], com_sat[:,2], vel_sat[:,0], vel_sat[:,1], vel_sat[:,2]]).T)
    np.savetxt(out_name_host, np.array([com_disk[:,0], com_disk[:,1], com_disk[:,2], vel_disk[:,0], vel_disk[:,1], vel_disk[:,2]]).T)



    
