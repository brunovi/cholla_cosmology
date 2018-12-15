import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import gc

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

def load_snapshot_data_grid( nSnap, inDir):
  inFileName = 'grid_{0}.h5'.format(nSnap)
  snapFile = h5.File( inDir + inFileName, 'r')
  inputKeys = snapFile.keys()
  grid_keys = [ 'density', 'momentum_x', 'momentum_y', 'momentum_z', 'GasEnergy']
  optional_keys = [ ]
  data_grid = {}
  for key in optional_keys:
    if key in inputKeys: grid_keys.append( key )
  for key in grid_keys:
    data_grid[key] = snapFile[key]
  return data_grid


inDir = '/lustre/atlas/proj-shared/ast125/data_cosmo/particles/'
inDir = '/lustre/atlas/proj-shared/ast125/data_cosmo/grid/'
outDir = inDir + 'projections/'
out_base_name = 'proj_ge'

nSnap = 0
for nSnap in range(20,30):
  data = load_snapshot_data_grid( nSnap, inDir )
  dens = data['density'][...]
  ge = data['GasEnergy'][...]
  dens2_proj = (ge*dens).sum( axis=0 ) / (dens).sum(axis=0)

  fileName = out_base_name + '_{0}.h5'.format( nSnap )
  fileSnap = h5.File( outDir + fileName, 'w' )

  fileSnap.create_dataset( 'density', data=dens2_proj )
  fileSnap.close()
  gc.collect()
  # data_particles = load_snapshot_data_particles( nSnap, inDir )
  # current_z = data_particles['current_z']
  # dens = data_particles['density'][...]
  # dens2_proj = (dens*dens*dens).sum( axis=0 ) / (dens*dens).sum(axis=0)
  #
  # fileName = out_base_name + '_{0}.h5'.format( nSnap )
  # fileSnap = h5.File( outDir + fileName, 'w' )
  #
  # fileSnap.attrs['current_z'] = current_z
  # fileSnap.create_dataset( 'density', data=dens2_proj )
  # fileSnap.close()

# box_size = 50
# d_min = dens2_proj.min()
# d_min = 0
# d_max = dens2_proj.max()

# plt.figure(0, figsize=(15,15))
# plt.clf()
# fig = plt.gcf()
# ax = plt.gca()
# img = np.log10( dens2_proj  )
# ax.imshow(img, extent=[0,box_size,0,box_size], interpolation='bilinear', cmap='inferno'  )
# ax.set_axis_off()
# text = r'z = {0:.1f}'.format(current_z)
# ax.text(box_size*0.05, box_size*0.90, text, fontsize=20, color='white',
#         bbox={'facecolor':'black', 'edgecolor':'white', 'alpha':0.5, })
#
# fig.axes[0].get_yaxis().set_visible(False)
# fig.axes[0].get_xaxis().set_visible(False)
# # # plt.axis('off')
# outFileName = outDir + 'densProj_{0}.png'.format(snap)
# fig.savefig( outFileName, pad_inches=0,  bbox_inches='tight', dpi=200 )
# # fig.show()
