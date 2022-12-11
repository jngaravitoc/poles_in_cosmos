"""
Script to read snapshots from cosmological simulations

@author:
  Nico Garavito-Camargo (github/jngaravitoc)
"""

import numpy as np
import gizmo_analysis as ga
import utilities as ut
import pynbody

def pynbody_snapshot(part_properties, part_data):
  """
  Families (partypes) in pynbody are defined here: https://github.com/pynbody/pynbody/blob/master/pynbody/default_config.ini

  Here we use will follow this convention:
    dm : all dark matter particles in halo
    dm_tracer : substructure in the dark matter halo
    star : all star particles in halo
    star_tracer : for stellar substructure 
    cloud : disk particles 
    
  """
  # get particle types
  
  # get particle lengths


  nfamilies = np.zeros(len(part_properties.keys()))
  i=0
  for particles in part_properties.keys():
    nfamilies[i] = length_family[particles+"_length"]
    family_name[i] = particles
    i+=1
  

  #h = pynbody.new(dark=, star=, debris=, cloud=, dm_tracer=, star_tracer=,  debris_tracer= )
  h=1
  """
  for particles in part_properties.keys():
    for properties in part_properties[particles]:
      print(particles, properties)
      pp = part[particles].prop(properties)
      new_property = {part+"_"+properties : pp}
      part_dict.update(new_propery)
  """
  return h 

def FIRE_snapshot(directory, snapname, particles_prop, **kwargs):
  """
  Snapshots species, units and more details can be found at:
  https://bitbucket.org/ngaravito/gizmo_analysis/src/master/gizmo_io.py
  
  TODO: Check units

  """

  part = ga.io.Read.read_snapshots(species=particles_prop['species'],
                                   snapshot_values=particles_prop['snapshot_values'],
                                   simulation_directory=particles_prop['simulation_directory'],
                                   simulation_name=particles_prop['simulation_name'],
                                   properties=particles_prop['properties'],
                                   particle_subsample_factor=particles_prop['particle_subsample_factor'],
                                   sort_dark_by_id =  particles_prop['sort_dark_by_id'],
                                   host_number = particles_prop['host_number'],
                                   assign_hosts = particles_prop['assign_hosts'],
                                   assign_orbits = particles_prop['assign_orbits'],
                                   assign_formation_coordinates = particles_prop['assign_formation_coordinates'])

  # Initialize particles_dictionary

  part_data = { 
      }
  species = particles_prop['species']
  properties = particles_prop['properties']
  
  for sp in species:
    for p in properties:
        print(sp, p)
        pp = part[sp].prop(p)
        new_prop = {sp+"_"+p : pp}
        part_data.update(new_prop)
  
  # do intermediate selections. Substructure, disk stars, etc..
  """
  if "disk" in particles:
    form_dist = part["star"].prop('form.host.distance.total')
    #form_time = part["star"].prop('form.time')
    disk_particles = np.where(form_dist<30)[0]
    pos_disk = stars_host_distance[disk_particles]
    vel_disk = stars_host_velocity[disk_particles]

  elif "dm_substructure" in particles:

  elif "stellar_substructure" in particles:
  """

  # pynbody halo 

  pynbody_dict = {
      }

  if "dark" in particles:
    npart = len(part[particles]["mass"])
    family_name = {"dark" : npart}
    pynbody_dict.update(family_name)

  elif "star" in particles:
    npart = len(part[particles]["mass"])
    family_name = {"star" : npart} 
    pynbody_dict.update(family_name)

  elif "disk" in particles:
    npart = len(part[particles]["mass"])
    family_name = {"cloud" : npart}
    pynbody_dict.update(family_name)

  elif "dm_substructure" in particles:
    npart = len(part[particles]["mass"])
    family_name = {"dm_tracer" : npart}
    pynbody_dict.update(family_name)

  elif "stellar_substructure" in particles:
    npart = len(part[particles]["mass"])
    family_name = {"star_tracer" : npart}
    pynbody_dict.update(family_name)

  else :
    print("Particle type not found")

  # make pynbody halo
  halo = pynbody_snapshot(pynbody_dict, part_data)

  return halo

def Auriga_snapshot(directory, snapname, partypes):
  return 0

def read_snapshot(simname, directory, snapnum, partypes, **kwargs):
  """
  simname : str 
    At the moment only FIRE or Auriga are supported
  partypes: ['star', 'dark', 'disk']
    {star: "positions", "masses", "velocities"}
    units? 

  """
  if simname == 'FIRE':
    snapshot = FIRE_snapshot(directory, snapnum, partypes, **kwargs)
  elif simname == "Auriga":
    snapshot = Auriga_snapshot(directory, snapnum, partypes)

  return snapshot



