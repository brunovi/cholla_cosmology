#ifdef PARTICLES
#include <unistd.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>
#ifdef HDF5
#include<hdf5.h>
#endif
#ifdef MPI_CHOLLA
#include"mpi_routines.h"
#endif

#include"io_particles.h"


// /* Write the initial conditions */
void WriteData_Particles( Grid3D &G, struct parameters P, int nfile)
{
// #ifndef  MPI_CHOLLA
  /*just use the data output routine*/
  OutputData_Particles( G, P, nfile);
// #else  /*MPI_CHOLLA*/
//   OutputDataMPI_Particles( G, P, nfile);
// #endif /*MPI_CHOLLA*/
}



#ifdef HDF5
/*! \fn void Write_Header_HDF5(hid_t file_id)
 *  \brief Write the relevant header info to the HDF5 file. */
void Write_Particles_Header_HDF5(Particles_3D &Particles, hid_t file_id)
{
  hid_t     attribute_id, dataspace_id;
  herr_t    status;
  hsize_t   attr_dims;
  int       int_data[3];
  Real      Real_data[3];

  // Single attributes first
  attr_dims = 1;
  // Create the data space for the attribute
  dataspace_id = H5Screate_simple(1, &attr_dims, NULL);
  // Create a group attribute
  attribute_id = H5Acreate(file_id, "t_particles", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // Write the attribute data
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Particles.t);
  // Close the attribute
  status = H5Aclose(attribute_id);
  attribute_id = H5Acreate(file_id, "dt_particles", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Particles.dt);
  status = H5Aclose(attribute_id);
  attribute_id = H5Acreate(file_id, "n_particles_local", H5T_STD_I64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_ULONG, &Particles.n_local);
  status = H5Aclose(attribute_id);
// #ifdef COSMOLOGY
//   attribute_id = H5Acreate(file_id, "current_a", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &current_a);
//   status = H5Aclose(attribute_id);
//   attribute_id = H5Acreate(file_id, "current_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &current_z);
//   status = H5Aclose(attribute_id);
// #endif
  status = H5Sclose(dataspace_id);
}

void Write_Particles_Data_HDF5( Particles_3D &Particles, hid_t file_id)
{
  part_int_t i, j, k, id, buf_id;
  hid_t     dataset_id, dataspace_id;
  Real      *dataset_buffer;
  part_int_t  *dataset_buffer_IDs;
  herr_t    status;
  part_int_t n_local = Particles.n_local;
  // int       nx_dset = H.nx_real;
  hsize_t   dims[1];
  dataset_buffer = (Real *) malloc(n_local*sizeof(Real));
  dataset_buffer_IDs = (part_int_t *) malloc(n_local*sizeof(part_int_t));

  // Create the data space for the datasets
  dims[0] = n_local;
  dataspace_id = H5Screate_simple(1, dims, NULL);


  // Copy the pos_x vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.pos_x[i];
  dataset_id = H5Dcreate(file_id, "/pos_x", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the pos_y vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.pos_y[i];
  dataset_id = H5Dcreate(file_id, "/pos_y", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the pos_z vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.pos_z[i];
  dataset_id = H5Dcreate(file_id, "/pos_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the vel_x vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.vel_x[i];
  dataset_id = H5Dcreate(file_id, "/vel_x", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the vel_y vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.vel_y[i];
  dataset_id = H5Dcreate(file_id, "/vel_y", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the vel_z vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.vel_z[i];
  dataset_id = H5Dcreate(file_id, "/vel_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the mass vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.mass[i];
  dataset_id = H5Dcreate(file_id, "/mass", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // 3D case
  int       nx_dset = Particles.G.nx_local;
  int       ny_dset = Particles.G.ny_local;
  int       nz_dset = Particles.G.nz_local;
  hsize_t   dims3d[3];
  dataset_buffer = (Real *) malloc(Particles.G.nz_local*Particles.G.ny_local*Particles.G.nx_local*sizeof(Real));

  // Create the data space for the datasets
  dims3d[0] = nx_dset;
  dims3d[1] = ny_dset;
  dims3d[2] = nz_dset;
  dataspace_id = H5Screate_simple(3, dims3d, NULL);

  // Copy the density array to the memory buffer
  int nGHST = Particles.G.n_ghost_particles_grid;
  for (k=0; k<Particles.G.nz_local; k++) {
    for (j=0; j<Particles.G.ny_local; j++) {
      for (i=0; i<Particles.G.nx_local; i++) {
        id = (i+nGHST) + (j+nGHST)*(Particles.G.nx_local+2*nGHST) + (k+nGHST)*(Particles.G.nx_local+2*nGHST)*(Particles.G.ny_local+2*nGHST);
        buf_id = k + j*Particles.G.nz_local + i*Particles.G.nz_local*Particles.G.ny_local;
        dataset_buffer[buf_id] = Particles.G.density[id];
      }
    }
  }

  // Create a dataset id for density
  dataset_id = H5Dcreate(file_id, "/density", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);

  free(dataset_buffer);
  free(dataset_buffer_IDs);
}


#endif

void OutputData_Particles( Grid3D &G, struct parameters P, int nfile)
{
  FILE *out;
  char filename[80];
  char timestep[20];

  // create the filename
  strcpy(filename, P.outdir);
  sprintf(timestep, "%d", nfile);
  strcat(filename,timestep);
  // a binary file is created for each process
  #if defined BINARY
  chprintf("\nERROR: Particles only support HDF5 outputs\n")
  return;
  // only one HDF5 file is created
  #elif defined HDF5
  strcat(filename,"_particles");
  strcat(filename,".h5");
  #ifdef MPI_CHOLLA
  sprintf(filename,"%s.%d",filename,procID);
  #endif
  #endif

  #if defined HDF5
  hid_t   file_id;
  herr_t  status;

  // Create a new file collectively
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Write header (file attributes)
  G.Write_Header_HDF5(file_id);
  Write_Particles_Header_HDF5( G.Particles, file_id);
  Write_Particles_Data_HDF5( G.Particles, file_id);

  // Close the file
  status = H5Fclose(file_id);
  #endif
}
#endif
