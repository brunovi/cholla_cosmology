import os, sys
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm


inDir = '/lustre/atlas/scratch/bvilasen/ast125/cosmo_1024_hydro/data/'
outDir = inDir + 'grid/'

name_base = 'h5'
out_base_name = 'grid_'

dataFiles = [f for f in listdir(inDir) if (isfile(join(inDir, f)) and (f.find('.h5.') > 0 ) and ( f.find('_particles') < 0) ) ]
dataFiles = np.sort( dataFiles )
nFiles = len( dataFiles )


def split_name( file_name):
  nSap, name, nBox = file_name.split('.')
  return [int(nSap), int(nBox)]

files_names = np.array([ split_name( file_name ) for file_name in dataFiles ])
snaps, boxes = files_names.T
snaps = np.unique( snaps )
boxes = np.unique( boxes )
nSnapshots = len( snaps )
nBoxes = len( boxes )

print "Number of boxes: {0}".format(nBoxes)
print "Number of snapshots: {0}".format(nSnapshots)

snap_id = 0
box_id = 0
inFileName = '{0}.{1}.{2}'.format(snap_id, name_base, box_id)
inFile = h5py.File( inDir + inFileName, 'r')
print inFile.keys()
head = inFile.attrs
dims_all = head['dims']
dims_local = head['dims_local']

keys = [  'density', 'momentum_x', 'momentum_y', 'momentum_z'  ]
if 'GasEnergy' in inFile.keys(): keys.append('GasEnergy')

nz, ny, nx = dims_all
inFile.close()

for nSnap in range(nSnapshots):
  fileName = out_base_name + '{0}.h5'.format( nSnap )
  fileSnap = h5py.File( outDir + fileName, 'w' )
  added_time = False
  print ' snap: {0}  {1}'.format( nSnap, keys )
  for key in keys:
    print '  {0}'.format(key)
    data_all = np.zeros( dims_all, dtype=np.float64 )
    for nBox in range(nBoxes):
      inFileName = '{0}.{1}.{2}'.format(nSnap, name_base, nBox)
      inFile = h5py.File( inDir + inFileName, 'r')
      head = inFile.attrs
      time = head['t']
      procStart_z, procStart_y, procStart_x = head['offset']
      procEnd_z, procEnd_y, procEnd_x = head['offset'] + head['dims_local']
      data_local = inFile[key][...]
      data_all[ procStart_z:procEnd_z, procStart_y:procEnd_y, procStart_x:procEnd_x] = data_local
      inFile.close()
    if key=='grav_density': print '  {0}   {1}   {2}'.format(data_all.mean(), data_all.min(), data_all.max())
    fileSnap.create_dataset( key, data=data_all.astype(np.float32) )
    fileSnap.attrs['max_'+ key ] = data_all.max()
    fileSnap.attrs['min_'+ key ] = data_all.min()
    fileSnap.attrs['mean_'+ key ] = data_all.mean()
    if added_time == False:
      fileSnap['t'] = time
      added_time = True
  fileSnap.close()
  print ' Saved File: ', outDir+fileName
