import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np


def load_snapshot_data_particles( nSnap, inputDir ):
  inFileName = 'particles_{0}.h5'.format(nSnap)
  partsFile = h5.File( inputDir + inFileName, 'r')
  fields_data = partsFile.keys()
  current_a = partsFile.attrs['current_a']
  current_z = partsFile.attrs['current_z']
  particle_mass = partsFile.attrs['particle_mass']

  data_part = {}
  data_part['current_a'] = current_a
  data_part['current_z'] = current_z
  data_part['particle_mass'] = particle_mass
  part_keys = [ 'density' ]
  extra_keys = [  ]
  for key in extra_keys:
    if key not in fields_data: continue
    if key in partsFile.keys(): part_keys.append(key)
  for key in part_keys:
    if key not in fields_data: continue
    data_part[key] = partsFile[key]
  return data_part

def load_snapshot_data_grid( nSnap, inFileName ):
  # inFileName = '{0}.h5'.format(nSnap)
  snapFile = h5.File( inFileName, 'r')
  t = snapFile.attrs['t'][0]
  inputKeys = snapFile.keys()
  grid_keys = [ 'density', 'momentum_x', 'momentum_y', 'momentum_z', 'GasEnergy']
  optional_keys = [ ]
  data_grid = {}
  data_grid['t'] = t
  for key in optional_keys:
    if key in inputKeys: grid_keys.append( key )
  for key in grid_keys:
    data_grid[key] = snapFile[key]
  return data_grid


inDir = '/lustre/atlas/proj-shared/ast125/data_cosmo/particles/'
outDir = inDir + 'projections/'

nSnap = 0


data_particles = load_snapshot_data_particles( nSnap, inDir )
