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
#include"io_particles.h"
#include"grid3D.h"
#ifdef MPI_CHOLLA
#include"mpi_routines.h"
#endif

#ifdef COSMOLOGY
#include "cosmology.h"
#endif

#ifdef TILING
#include "error_handling.h"
#endif //TILING

// /* Write the initial conditions */
void WriteData_Particles(Part3D &Particles, Grid3D G, struct parameters P, int nfile)
{
#ifndef  MPI_CHOLLA
  /*just use the data output routine*/
  OutputData_Particles( Particles, G,P,nfile);
#else  /*MPI_CHOLLA*/
  OutputDataMPI_Particles( Particles, G,P,nfile);
#endif /*MPI_CHOLLA*/
}


#ifdef MPI_CHOLLA
/* Output the grid data to file. */
void OutputDataMPI_Particles(Part3D &Particles, Grid3D G, struct parameters P, int nfile)
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
  strcat(filename,".bin");
  sprintf(filename,"%s.%d",filename,procID);
  // only one HDF5 file is created
  #elif defined HDF5
  strcat(filename,"_parts");
  strcat(filename,".h5");
  sprintf(filename,"%s.%d",filename,procID);
  #endif

  // open the files for binary writes
  #if defined BINARY
  out = fopen(filename, "w");
  if(out == NULL) {printf("Error opening output file.\n"); fflush(stdout); exit(0); }

  // write the header to the output files
  G.Write_Header_Binary(out);

  // write the conserved variables to the output files
  G.Write_Grid_Binary(out);

  // close the output file
  fclose(out);


  #elif defined HDF5
  hid_t   file_id;
  herr_t  status;

  // Create a new file collectively
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Write header (file attributes)
  G.Write_Header_HDF5(file_id);

  // Write the conserved variables to the output file
  // G.Write_Grid_HDF5(file_id);

  #ifdef PARTICLES
  Particles.Write_PartHeader_HDF5(file_id);
  Particles.Write_PartData_HDF5(file_id);
  #endif
  // Close the file
  status = H5Fclose(file_id);
  #endif
}
#endif /*MPI_CHOLLA*/




/*! \fn void Load_PartData(struct parameters P)
 *  \brief Read in particle data from an output file. */
void Part3D::Load_PartData(struct parameters P) {

  char filename[100];
  // char timestep[20];
  // int nfile = P.nfile; //output step you want to read from
  char filename_counter[100];
  // create the filename to read from
  if ( P.nfile == 0){
    // strcpy(filename, P.indir);
    // // assumes your data is in the outdir specified in the input file
    // // strcat(filename,"/home/bruno/Desktop/Dropbox/Developer/cholla_analysis/data/cosmo_ics/ics_256_f");
    // // strcat(filename,"/home/bruno/Desktop/Dropbox/Developer/cholla_analysis/data/cosmo_ics/cosmo_ICs_128_float");
    // // strcat(filename,"/raid/bruno/data/cosmo_sims/cholla_pm/ics/cosmo_ICs_256_float");
    // strcat(filename,"snapshot_000_particles");

    strcpy(filename, P.indir);
    strcat(filename,"0_parts");
  }
  else{
    char timestep[20];
    int nfile = P.nfile; //output step you want to read from

    // nfile = 0; //BRUNO WAS HERE

    // create the filename to read from
    // assumes your data is in the outdir specified in the input file
    // strcpy(filename, P.outdir);
    strcpy(filename, P.indir);
    sprintf(timestep, "%d_parts", nfile);
    strcat(filename,timestep);
  }


  // strcpy(filename, P.outdir);
  // sprintf(timestep, "%d", nfile);
  // strcat(filename,timestep);
  #if defined BINARY
  strcat(filename,".bin");
  #elif defined HDF5
  strcat(filename,".h5");
  #endif
  // for now assumes you will run on the same number of processors
  #ifdef MPI_CHOLLA
  if ( P.nfile > 0) sprintf(filename,"%s.%d",filename,procID);
  #endif

  if ( P.nfile > 0 ){
    chprintf("Loading Particle data: %s\n", filename);

    #if defined BINARY
    FILE *fp;
    // open the file
    fp = fopen(filename, "r");
    if (!fp) {
      printf("Unable to open input file.\n");
      exit(0);
    }

    // read in grid data
    Read_Grid_Binary(fp);

    // close the file
    fclose(fp);

    #elif defined HDF5
    hid_t  file_id;
    herr_t  status;

    // open the file
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
      printf("Unable to open input file.\n");
      exit(0);
    }

    // read in grid data
    Load_PartData_HDF5(file_id, P.nfile, true);

    // close the file
    status = H5Fclose(file_id);
    #endif
  }

  else{
    int counter, file_number;
    int n_files_particles = P.n_parts_initFiles;
    chprintf( "\nLoading particles from %d files...\n", n_files_particles);
    for( counter = 0; counter<n_files_particles; counter++ ){
      // file_number = counter + procID;
      // if ( file_number >= nproc ) file_number -= nproc;
      file_number = counter; 

      sprintf(filename_counter,"%s.%d",filename,file_number);
      chprintf( "Loading File: %s\n", filename_counter);

      hid_t  file_id;
      herr_t  status;

      // open the file
      file_id = H5Fopen(filename_counter, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file_id < 0) {
        printf("Unable to open input file.\n");
        exit(0);
      }

      // read in grid data
      if (counter == n_files_particles-1 ) Load_PartData_HDF5(file_id, P.nfile, true);
      else Load_PartData_HDF5(file_id, P.nfile, false);

      // close the file
      status = H5Fclose(file_id);
    }
  }

}









#ifdef HDF5
/*! \fn void Write_Header_HDF5(hid_t file_id)
 *  \brief Write the relevant header info to the HDF5 file. */
void Part3D::Write_PartHeader_HDF5(hid_t file_id)
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
  attribute_id = H5Acreate(file_id, "t_parts", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // Write the attribute data
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &t);
  // Close the attribute
  status = H5Aclose(attribute_id);
  attribute_id = H5Acreate(file_id, "dt_parts", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &dt);
  status = H5Aclose(attribute_id);
  attribute_id = H5Acreate(file_id, "nParts_local", H5T_STD_I64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_ULONG, &n_local);
  status = H5Aclose(attribute_id);
#ifdef COSMOLOGY
  attribute_id = H5Acreate(file_id, "current_a", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &current_a);
  status = H5Aclose(attribute_id);
  attribute_id = H5Acreate(file_id, "current_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &current_z);
  status = H5Aclose(attribute_id);
#endif
  // attribute_id = H5Acreate(file_id, "nParts_total", H5T_STD_I64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // status = H5Awrite(attribute_id, H5T_NATIVE_ULONG, &n_local);
  // status = H5Aclose(attribute_id);
  // Close the dataspace
  status = H5Sclose(dataspace_id);

  // // Now 3D attributes
  // attr_dims = 3;
  // // Create the data space for the attribute
  // dataspace_id = H5Screate_simple(1, &attr_dims, NULL);
  //
  // #ifndef MPI_CHOLLA
  // int_data[0] = H.nx_real;
  // int_data[1] = H.ny_real;
  // int_data[2] = H.nz_real;
  // #endif
  // #ifdef MPI_CHOLLA
  // int_data[0] = nx_global_real;
  // int_data[1] = ny_global_real;
  // int_data[2] = nz_global_real;
  // #endif
  //
  // attribute_id = H5Acreate(file_id, "dims", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // status = H5Awrite(attribute_id, H5T_NATIVE_INT, int_data);
  // status = H5Aclose(attribute_id);
  //
  // #ifdef MPI_CHOLLA
  // int_data[0] = H.nx_real;
  // int_data[1] = H.ny_real;
  // int_data[2] = H.nz_real;
  //
  // attribute_id = H5Acreate(file_id, "dims_local", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // status = H5Awrite(attribute_id, H5T_NATIVE_INT, int_data);
  // status = H5Aclose(attribute_id);
  //
  // int_data[0] = nx_local_start;
  // int_data[1] = ny_local_start;
  // int_data[2] = nz_local_start;
  //
  // attribute_id = H5Acreate(file_id, "offset", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // status = H5Awrite(attribute_id, H5T_NATIVE_INT, int_data);
  // status = H5Aclose(attribute_id);
  // #endif
  //
  // Real_data[0] = H.xbound;
  // Real_data[1] = H.ybound;
  // Real_data[2] = H.zbound;
  //
  // attribute_id = H5Acreate(file_id, "bounds", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, Real_data);
  // status = H5Aclose(attribute_id);
  //
  // Real_data[0] = H.xdglobal;
  // Real_data[1] = H.ydglobal;
  // Real_data[2] = H.zdglobal;
  //
  // attribute_id = H5Acreate(file_id, "domain", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, Real_data);
  // status = H5Aclose(attribute_id);
  //
  // Real_data[0] = H.dx;
  // Real_data[1] = H.dy;
  // Real_data[2] = H.dz;
  //
  // attribute_id = H5Acreate(file_id, "dx", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  // status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, Real_data);
  // status = H5Aclose(attribute_id);
  //
  // // Close the dataspace
  // status = H5Sclose(dataspace_id);

}
/*! \fn void Write_Grid_HDF5(hid_t file_id)
 *  \brief Write the grid to a file, at the current simulation time. */
void Part3D::Write_PartData_HDF5(hid_t file_id)
{
  partID_t i, j, k, id, buf_id;
  hid_t     dataset_id, dataspace_id;
  Real      *dataset_buffer;
  partID_t  *dataset_buffer_IDs;
  herr_t    status;

  // int       nx_dset = H.nx_real;
  hsize_t   dims[1];
  dataset_buffer = (Real *) malloc(n_local*sizeof(Real));
  dataset_buffer_IDs = (partID_t *) malloc(n_local*sizeof(partID_t));

  // Create the data space for the datasets
  dims[0] = n_local;
  dataspace_id = H5Screate_simple(1, dims, NULL);


  // // Copy the pos_x vector to the memory buffer
  // for ( i=0; i<n_local; i++) dataset_buffer_IDs[i] = partIDs[i];
  // // Create a dataset id for pos_x
  // dataset_id = H5Dcreate(file_id, "/partID", H5T_STD_I64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  // status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_IDs);
  // // Free the dataset id
  // status = H5Dclose(dataset_id);



  // Copy the pos_x vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = pos_x[i];
  // Create a dataset id for pos_x
  dataset_id = H5Dcreate(file_id, "/pos_x", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);

  // Copy the pos_y vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = pos_y[i];
  // Create a dataset id for pos_x
  dataset_id = H5Dcreate(file_id, "/pos_y", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);

  // Copy the pos_z vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = pos_z[i];
  // Create a dataset id for pos_x
  dataset_id = H5Dcreate(file_id, "/pos_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);

  Real vel_norm = 1.0;
#ifdef COSMOLOGY
  // vel_norm = v_0 /( current_a * sqrt(current_a));  //From GADGET guide to get peculiar velocities
  vel_norm = 1;
#endif
  // Copy the vel_x vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = vel_x[i]*vel_norm;
  dataset_id = H5Dcreate(file_id, "/vel_x", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the vel_y vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = vel_y[i]*vel_norm;
  dataset_id = H5Dcreate(file_id, "/vel_y", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the vel_z vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = vel_z[i]*vel_norm;
  dataset_id = H5Dcreate(file_id, "/vel_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);

  // Copy the mass vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = mass[i];
  dataset_id = H5Dcreate(file_id, "/mass", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);
  //
  // status = H5Sclose(dataspace_id);

  // 3D case

  int       nx_dset = G.nx;
  int       ny_dset = G.ny;
  int       nz_dset = G.nz;
  hsize_t   dims3d[3];
  dataset_buffer = (Real *) malloc(G.nz*G.ny*G.nx*sizeof(Real));

  // Create the data space for the datasets
  dims3d[0] = nx_dset;
  dims3d[1] = ny_dset;
  dims3d[2] = nz_dset;
  dataspace_id = H5Screate_simple(3, dims3d, NULL);

  // Copy the density array to the memory buffer
  int nGHST = G.grid_ghost;
  for (k=0; k<G.nz; k++) {
    for (j=0; j<G.ny; j++) {
      for (i=0; i<G.nx; i++) {
        id = (i+nGHST) + (j+nGHST)*(G.nx+2*nGHST) + (k+nGHST)*(G.nx+2*nGHST)*(G.ny+2*nGHST);
        buf_id = k + j*G.nz + i*G.nz*G.ny;
        dataset_buffer[buf_id] = G.density_h[id];
      }
    }
  }

  // Create a dataset id for density
  dataset_id = H5Dcreate(file_id, "/density", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);




  // Copy the density array to the memory buffer
  int nGHST_pot = G.grid_ghost_pot;
  for (k=0; k<G.nz; k++) {
    for (j=0; j<G.ny; j++) {
      for (i=0; i<G.nx; i++) {
        id = (i+nGHST_pot) + (j+nGHST_pot)*(G.nx+2*nGHST_pot) + (k+nGHST_pot)*(G.nx+2*nGHST_pot)*(G.ny+2*nGHST_pot);
        buf_id = k + j*G.nz + i*G.nz*G.ny;
        dataset_buffer[buf_id] = G.potential[id];
      }
    }
  }

  // Create a dataset id for density
  dataset_id = H5Dcreate(file_id, "/potential", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);


  // // Free the dataspace id
  // status = H5Sclose(dataspace_id);


  free(dataset_buffer);
  free(dataset_buffer_IDs);

}


void Part3D::Load_PartData_HDF5(hid_t file_id, int nfile, bool print_load )
{
  int i, j, k, id, buf_id;
  hid_t     attribute_id, dataset_id;
  Real      *dataset_buffer_px;
  Real      *dataset_buffer_py;
  Real      *dataset_buffer_pz;
  Real      *dataset_buffer_vx;
  Real      *dataset_buffer_vy;
  Real      *dataset_buffer_vz;
  Real      *dataset_buffer_m;
  herr_t    status;

  partID_t N_DM, pIndx;
  partID_t N_DM_Total;
  //
  // Read in header values not set by grid initialization
  if ( nfile == 0 ){
    attribute_id = H5Aopen(file_id, "N_DM", H5P_DEFAULT);
    status = H5Aread(attribute_id, H5T_NATIVE_LONG, &N_DM_Total);
    status = H5Aclose(attribute_id);
#ifdef TILING
    N_DM = N_DM_Total;
    N_DM_Total *= nproc;
#else  //TILING
    attribute_id = H5Aopen(file_id, "N_DM_file", H5P_DEFAULT);
    status = H5Aread(attribute_id, H5T_NATIVE_LONG, &N_DM);
    status = H5Aclose(attribute_id);
#endif //TILING
    chprintf(" Loading %ld  total particles...\n", N_DM_Total );
    chprintf("  Reading %ld  in file particles...\n", N_DM );
  }
  else{
    attribute_id = H5Aopen(file_id, "nParts_local", H5P_DEFAULT);
    status = H5Aread(attribute_id, H5T_NATIVE_LONG, &N_DM);
    status = H5Aclose(attribute_id);
    attribute_id = H5Aopen(file_id, "current_a", H5P_DEFAULT);
    status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &current_a);
    status = H5Aclose(attribute_id);
    attribute_id = H5Aopen(file_id, "current_z", H5P_DEFAULT);
    status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &current_z);
    status = H5Aclose(attribute_id);
    chprintf(" Loading %ld particles...\n", N_DM );
    chprintf(" Current scale_factor: %.4f \n", current_a );
    current_a_gas = current_a;
  }

  // attribute_id = H5Aopen(file_id, "t", H5P_DEFAULT);
  // status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &H.t);
  // status = H5Aclose(attribute_id);
  // attribute_id = H5Aopen(file_id, "dt", H5P_DEFAULT);
  // status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &H.dt);
  // status = H5Aclose(attribute_id);
  // attribute_id = H5Aopen(file_id, "n_step", H5P_DEFAULT);
  // status = H5Aread(attribute_id, H5T_NATIVE_INT, &H.n_step);
  // status = H5Aclose(attribute_id);
  //
  // // 1D case
  // if (H.nx>1 && H.ny==1 && H.nz==1) {
  //
  // need a dataset buffer to remap fastest index
  //dataset_buffer_px = (Real *) malloc(N_DM*sizeof(Real));
  if(!(dataset_buffer_px = (Real *) malloc(N_DM*sizeof(Real))))
  {
    printf("Error allocating dataset_buffer_px on procID %d\n",procID);
  }

  // Open the pos_x dataset
  dataset_id = H5Dopen(file_id, "/pos_x", H5P_DEFAULT);
  // Read the density array into the dataset buffer // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_px);
  // Free the dataset id
  status = H5Dclose(dataset_id);

  chprintf("Read x\n");

  if(!(dataset_buffer_py = (Real *) malloc(N_DM*sizeof(Real))))
  {
    printf("Error allocating dataset_buffer_py on procID %d\n",procID);
  }

  dataset_id = H5Dopen(file_id, "/pos_y", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_py);
  status = H5Dclose(dataset_id);

  chprintf("Read y\n");
  if(!(dataset_buffer_pz = (Real *) malloc(N_DM*sizeof(Real))))
  {
    printf("Error allocating dataset_buffer_pz on procID %d\n",procID);
  }

  dataset_id = H5Dopen(file_id, "/pos_z", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_pz);
  status = H5Dclose(dataset_id);

  chprintf("Read z\n");
  if(!(dataset_buffer_vx = (Real *) malloc(N_DM*sizeof(Real))))
  {
    printf("Error allocating dataset_buffer_vx on procID %d\n",procID);
  }
  dataset_id = H5Dopen(file_id, "/vel_x", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_vx);
  status = H5Dclose(dataset_id);

  chprintf("Read vx\n");
  if(!(dataset_buffer_vy = (Real *) malloc(N_DM*sizeof(Real))))
  {
    printf("Error allocating dataset_buffer_vy on procID %d\n",procID);
  }
  dataset_id = H5Dopen(file_id, "/vel_y", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_vy);
  status = H5Dclose(dataset_id);

  chprintf("Read vy\n");
  if(!(dataset_buffer_vz = (Real *) malloc(N_DM*sizeof(Real))))
  {
    printf("Error allocating dataset_buffer_vz on procID %d\n",procID);
  }
  dataset_id = H5Dopen(file_id, "/vel_z", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_vz);
  status = H5Dclose(dataset_id);

  chprintf("Read vz\n");
  if(!(dataset_buffer_m = (Real *) malloc(N_DM*sizeof(Real))))
  {
    printf("Error allocating dataset_buffer_m on procID %d\n",procID);
  }
  dataset_id = H5Dopen(file_id, "/mass", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_m);
  status = H5Dclose(dataset_id);

  chprintf("Read m\n");
  Real pPos_x=0, pPos_y=0, pPos_z=0;
  Real pVel_x=0, pVel_y=0, pVel_z=0, pMass=0;
  Real pos_norm = 1.0;
  Real vel_norm = 1.0;
#ifdef COSMO_UNITS
  pos_norm = 1;
  // vel_norm = v_0 /( current_a * sqrt(current_a));  //From GADGET guide to get peculiar velocities
  vel_norm = 1;  //From GADGET guide to get peculiar velocities
  // if ( nfile > 0 ) vel_norm = 1/( current_a * sqrt(current_a));  //From GADGET guide to get peculiar velocities
#endif
  Real mass_norm = 1.0;
  bool in_local;
  //printf("procID %d xMin %e yMin %e zMin %e\n",procID,G.xMin,G.yMin,G.zMin);
  //chprintf("xMin %e yMin %e zMin %e\n",G.xMin,G.yMin,G.zMin);
  for( pIndx=0; pIndx<N_DM; pIndx++ ){

#ifdef TILING
    //ERROR SOMETIMES HERE
    //pPos_x = dataset_buffer_px[pIndx] + G.xMin;
    //pPos_y = dataset_buffer_py[pIndx] + G.yMin;
    //pPos_z = dataset_buffer_pz[pIndx] + G.zMin;
    pPos_x = dataset_buffer_px[pIndx];
    pPos_y = dataset_buffer_py[pIndx];
    pPos_z = dataset_buffer_pz[pIndx];
#else //TILING
    pPos_x = dataset_buffer_px[pIndx];
    pPos_y = dataset_buffer_py[pIndx];
    pPos_z = dataset_buffer_pz[pIndx];
#endif //TILING

    pPos_x = ( pPos_x  ) / pos_norm;
    pPos_y = ( pPos_y  ) / pos_norm;
    pPos_z = ( pPos_z  ) / pos_norm;
    pVel_x = dataset_buffer_vx[pIndx] / vel_norm;
    pVel_y = dataset_buffer_vy[pIndx] / vel_norm;
    pVel_z = dataset_buffer_vz[pIndx] / vel_norm;
    pMass = dataset_buffer_m[pIndx] / mass_norm;
    // chprintf("%f\n", pVel_x*pVel_x + pVel_y*pVel_y + pVel_z*pVel_z  );
    //
    //
    //
/*
    in_local = true;
    // chprintf("%.02f \n", pPos_x);
    if ( pPos_x < G.domainMin_x || pPos_x > G.domainMax_x ){
      std::cout << " Particle outside global domain " << std::endl;
    }
    if ( pPos_y < G.domainMin_y || pPos_y > G.domainMax_y ){
      std::cout << " Particle outside global domain " << std::endl;
    }
    if ( pPos_z < G.domainMin_z || pPos_z > G.domainMax_z ){
      std::cout << " Particle outside global domain " << std::endl;
    }
    if ( pPos_x < G.xMin || pPos_x >= G.xMax ) in_local = false;
    if ( pPos_y < G.yMin || pPos_y >= G.yMax ) in_local = false;
    if ( pPos_z < G.zMin || pPos_z >= G.zMax ) in_local = false;
    // if ( G.xMax == G.domainMax_x && pPos_x == G.domainMax_x ) in_local = true;
    // if ( G.yMax == G.domainMax_y && pPos_y == G.domainMax_y ) in_local = true;
    // if ( G.zMax == G.domainMax_z && pPos_z == G.domainMax_z ) in_local = true;
    if ( ! in_local ) continue;
*/

#ifdef TILING
    partIDs.push_back(pIndx + procID*N_DM); //prevent duplicate IDs
#else //TILING
    partIDs.push_back(pIndx);    
#endif //TILING
    mass.push_back( pMass );
   // pos_x.push_back( pPos_x*pos_norm );
    //pos_y.push_back( pPos_y*pos_norm );
    //pos_z.push_back( pPos_z*pos_norm );
    pos_x.push_back( pPos_x*pos_norm + G.xMin);
    pos_y.push_back( pPos_y*pos_norm + G.yMin);
    pos_z.push_back( pPos_z*pos_norm + G.zMin);
    vel_x.push_back( pVel_x );
    vel_y.push_back( pVel_y );
    vel_z.push_back( pVel_z );
    grav_x.push_back( 0.0 );
    grav_y.push_back( 0.0 );
    grav_z.push_back( 0.0 );
    n_local += 1;

  }

  fflush(stdout);
  MPI_Barrier(world);
  chprintf("Particles read.\n");
  //
  //if (print_load){
  /*
  if (1){
    for(int i=0;i<nproc;i++)
    {
      if(i==procID)
      {
        printf("  [procID %d] Loaded: %ld \n",procID, n_local );
      }
      fflush(stdout);
      usleep(10);
      MPI_Barrier(world);
    }
  }
  */
  free(dataset_buffer_px);
  free(dataset_buffer_py);
  free(dataset_buffer_pz);
  free(dataset_buffer_vx);
  free(dataset_buffer_vy);
  free(dataset_buffer_vz);
  free(dataset_buffer_m);


  //chexit(0);
}
















#endif //HDF5
#endif //PARTICLES
