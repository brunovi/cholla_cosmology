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
#include"../mpi_routines.h"
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
  #ifdef COSMOLOGY
  attribute_id = H5Acreate(file_id, "current_a", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Particles.current_a);
  status = H5Aclose(attribute_id);
  attribute_id = H5Acreate(file_id, "current_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Particles.current_z);
  status = H5Aclose(attribute_id);
  #endif

  #ifdef SINGLE_PARTICLE_MASS
  attribute_id = H5Acreate(file_id, "particle_mass", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Particles.particle_mass);
  status = H5Aclose(attribute_id);
  #endif

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

  #ifndef SINGLE_PARTICLE_MASS
  // Copy the mass vector to the memory buffer
  for ( i=0; i<n_local; i++) dataset_buffer[i] = Particles.mass[i];
  dataset_id = H5Dcreate(file_id, "/mass", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);
  #endif

  #ifdef PARTICLE_IDS
  dataset_buffer_IDs = (part_int_t *) malloc(n_local*sizeof(part_int_t));
  for ( i=0; i<n_local; i++) dataset_buffer_IDs[i] = Particles.partIDs[i];
  dataset_id = H5Dcreate(file_id, "/particle_IDs", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_IDs);
  status = H5Dclose(dataset_id);
  free(dataset_buffer_IDs);
  #endif

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


void Load_Particles_Data( Particles_3D &Particles, struct parameters P){
  char filename[100];
  char timestep[20];
  int nfile = P.nfile; //output step you want to read from
  char filename_counter[100];
  // create the filename to read from

  strcpy(filename, P.indir);
  sprintf(timestep, "%d_particles", nfile);
  strcat(filename,timestep);

  #if defined BINARY
  chprintf("\nERROR: Particles only support HDF5 outputs\n")
  #elif defined HDF5
  strcat(filename,".h5");
  #endif

  #ifdef MPI_CHOLLA
  sprintf(filename,"%s.%d",filename,procID);
  #endif

  chprintf(" Loading particles file: %s \n", filename );

  #ifdef HDF5
  hid_t  file_id;
  herr_t  status;

  // open the file
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    printf("Unable to open input file.\n");
    exit(0);
  }

  Load_Particles_Data_HDF5(file_id, nfile, Particles );

  #endif
}

#ifdef HDF5
void Load_Particles_Data_HDF5(hid_t file_id, int nfile, Particles_3D &Particles )
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

  part_int_t n_local, pIndx;
  // Real part_dt, part_t;

  attribute_id = H5Aopen(file_id, "n_particles_local", H5P_DEFAULT);
  status = H5Aread(attribute_id, H5T_NATIVE_LONG, &n_local);
  status = H5Aclose(attribute_id);

  #ifdef COSMOLOGY
  attribute_id = H5Aopen(file_id, "current_z", H5P_DEFAULT);
  status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Particles.current_z);
  status = H5Aclose(attribute_id);

  attribute_id = H5Aopen(file_id, "current_a", H5P_DEFAULT);
  status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Particles.current_a);
  status = H5Aclose(attribute_id);
  #endif

  #ifdef SINGLE_PARTICLE_MASS
  attribute_id = H5Aopen(file_id, "particle_mass", H5P_DEFAULT);
  status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Particles.particle_mass);
  status = H5Aclose(attribute_id);
  chprintf( " Using Single mass for DM particles: %f  Msun/h\n", Particles.particle_mass);
  #endif

  #ifndef MPI_CHOLLA
  // if ( n_total != G.Grav.nx_total * G.Grav.ny_total * G.Grav.nz_total) break;
  chprintf(" Loading %ld particles\n", n_local);
  // #endif
  #else
  part_int_t n_total_load;
  n_total_load = ReducePartIntSum( n_local );
  chprintf( " Total Particles To Load: %ld\n", n_total_load );
  // for ( int i=0; i<nproc; i++ ){
  //   if ( procID == i ) std::cout << "  [pId:"  << procID << "]  Loading Particles: " << n_local <<  std::endl;
  //   MPI_Barrier(world);
  // }
  MPI_Barrier(world);
  #endif


  dataset_buffer_px = (Real *) malloc(n_local*sizeof(Real));
  dataset_id = H5Dopen(file_id, "/pos_x", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_px);
  status = H5Dclose(dataset_id);

  dataset_buffer_py = (Real *) malloc(n_local*sizeof(Real));
  dataset_id = H5Dopen(file_id, "/pos_y", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_py);
  status = H5Dclose(dataset_id);

  dataset_buffer_pz = (Real *) malloc(n_local*sizeof(Real));
  dataset_id = H5Dopen(file_id, "/pos_z", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_pz);
  status = H5Dclose(dataset_id);

  dataset_buffer_vx = (Real *) malloc(n_local*sizeof(Real));
  dataset_id = H5Dopen(file_id, "/vel_x", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_vx);
  status = H5Dclose(dataset_id);

  dataset_buffer_vy = (Real *) malloc(n_local*sizeof(Real));
  dataset_id = H5Dopen(file_id, "/vel_y", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_vy);
  status = H5Dclose(dataset_id);

  dataset_buffer_vz = (Real *) malloc(n_local*sizeof(Real));
  dataset_id = H5Dopen(file_id, "/vel_z", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_vz);
  status = H5Dclose(dataset_id);

  #ifndef SINGLE_PARTICLE_MASS
  dataset_buffer_m = (Real *) malloc(n_local*sizeof(Real));
  dataset_id = H5Dopen(file_id, "/mass", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_m);
  status = H5Dclose(dataset_id);
  #endif

  #ifdef PARTICLE_IDS
  part_int_t *dataset_buffer_IDs;
  dataset_buffer_IDs = (part_int_t *) malloc(n_local*sizeof(part_int_t));
  dataset_id = H5Dopen(file_id, "/particle_IDs", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer_IDs);
  status = H5Dclose(dataset_id);
  #endif

  Real px_min, px_max;
  Real py_min, py_max;
  Real pz_min, pz_max;
  Real vx_min, vx_max;
  Real vy_min, vy_max;
  Real vz_min, vz_max;
  px_min = 1e64;
  py_min = 1e64;
  pz_min = 1e64;
  px_max = -1e64;
  py_max = -1e64;
  pz_max = -1e64;



  Real pPos_x, pPos_y, pPos_z;
  Real pVel_x, pVel_y, pVel_z, pMass;
  part_int_t pID;
  bool in_local;
  for( pIndx=0; pIndx<n_local; pIndx++ ){
    pPos_x = dataset_buffer_px[pIndx];
    pPos_y = dataset_buffer_py[pIndx];
    pPos_z = dataset_buffer_pz[pIndx];
    pVel_x = dataset_buffer_vx[pIndx];
    pVel_y = dataset_buffer_vy[pIndx];
    pVel_z = dataset_buffer_vz[pIndx];
    #ifndef SINGLE_PARTICLE_MASS
    pMass = dataset_buffer_m[pIndx];
    #endif
    #ifdef PARTICLE_IDS
    pID = dataset_buffer_IDs[pIndx];
    #endif
    in_local = true;
    if ( pPos_x < Particles.G.domainMin_x || pPos_x > Particles.G.domainMax_x ){
      std::cout << " Particle outside global domain " << std::endl;
    }
    if ( pPos_y < Particles.G.domainMin_y || pPos_y > Particles.G.domainMax_y ){
      std::cout << " Particle outside global domain " << std::endl;
    }
    if ( pPos_z < Particles.G.domainMin_z || pPos_z > Particles.G.domainMax_z ){
      std::cout << " Particle outside global domain " << std::endl;
    }
    if ( pPos_x < Particles.G.xMin || pPos_x >= Particles.G.xMax ) in_local = false;
    if ( pPos_y < Particles.G.yMin || pPos_y >= Particles.G.yMax ) in_local = false;
    if ( pPos_z < Particles.G.zMin || pPos_z >= Particles.G.zMax ) in_local = false;
    if ( ! in_local ) {
      #ifdef PARTICLE_IDS
      std::cout << " Particle outside Loacal  domain    pID: " << pID << std::endl;
      #else
      std::cout << " Particle outside Loacal  domain " << std::endl;
      #endif
      std::cout << "  Domain X: " << Particles.G.xMin <<  "  " << Particles.G.xMax << std::endl;
      std::cout << "  Domain Y: " << Particles.G.yMin <<  "  " << Particles.G.yMax << std::endl;
      std::cout << "  Domain Z: " << Particles.G.zMin <<  "  " << Particles.G.zMax << std::endl;
      std::cout << "  Particle X: " << pPos_x << std::endl;
      std::cout << "  Particle Y: " << pPos_y << std::endl;
      std::cout << "  Particle Z: " << pPos_z << std::endl;
      continue;
    }

    if  ( pPos_x > px_max ) px_max = pPos_x;
    if  ( pPos_y > py_max ) py_max = pPos_y;
    if  ( pPos_z > pz_max ) pz_max = pPos_z;

    if  ( pPos_x < px_min ) px_min = pPos_x;
    if  ( pPos_y < py_min ) py_min = pPos_y;
    if  ( pPos_z < pz_min ) pz_min = pPos_z;

    if  ( pVel_x > vx_max ) vx_max = pVel_x;
    if  ( pVel_y > vy_max ) vy_max = pVel_y;
    if  ( pVel_z > vz_max ) vz_max = pVel_z;

    if  ( pVel_x < vx_min ) vx_min = pVel_x;
    if  ( pVel_y < vy_min ) vy_min = pVel_y;
    if  ( pVel_z < vz_min ) vz_min = pVel_z;

    Particles.pos_x.push_back( pPos_x );
    Particles.pos_y.push_back( pPos_y );
    Particles.pos_z.push_back( pPos_z );
    Particles.vel_x.push_back( pVel_x );
    Particles.vel_y.push_back( pVel_y );
    Particles.vel_z.push_back( pVel_z );
    Particles.grav_x.push_back( 0.0 );
    Particles.grav_y.push_back( 0.0 );
    Particles.grav_z.push_back( 0.0 );
    #ifndef SINGLE_PARTICLE_MASS
    Particles.mass.push_back( pMass );
    #endif
    #ifdef PARTICLE_IDS
    Particles.partIDs.push_back(pID);
    #endif
    Particles.n_local += 1;
  }
  #ifndef MPI_CHOLLA
  chprintf( " Loaded  %ld  particles\n", Particles.n_local );
  #else
  // for ( int i=0; i<nproc; i++ ){
  //   if ( procID == i ) std::cout << "  [pId:"  << procID << "]  N Particles Loaded: " << Particles.n_local <<  std::endl;
  //   MPI_Barrier(world);
  // }
  MPI_Barrier(world);
  part_int_t n_total_loaded;
  n_total_loaded = ReducePartIntSum( Particles.n_local );
  chprintf( " Total Particles Loaded: %ld\n", n_total_loaded );
  #endif

  #ifdef MPI_CHOLLA
  Real px_max_g = ReduceRealMax( px_max );
  Real py_max_g = ReduceRealMax( py_max );
  Real pz_max_g = ReduceRealMax( pz_max );
  Real vx_max_g = ReduceRealMax( vx_max );
  Real vy_max_g = ReduceRealMax( vy_max );
  Real vz_max_g = ReduceRealMax( vz_max );

  Real px_min_g = ReduceRealMax( px_min );
  Real py_min_g = ReduceRealMin( py_min );
  Real pz_min_g = ReduceRealMin( pz_min );
  Real vx_min_g = ReduceRealMin( vx_min );
  Real vy_min_g = ReduceRealMin( vy_min );
  Real vz_min_g = ReduceRealMin( vz_min );

  chprintf( "  Pos X   Min: %f   Max: %f   [ kpc/h]\n", px_min_g, px_max_g);
  chprintf( "  Pos Y   Min: %f   Max: %f   [ kpc/h]\n", py_min_g, py_max_g);
  chprintf( "  Pos Z   Min: %f   Max: %f   [ kpc/h]\n", pz_min_g, pz_max_g);
  chprintf( "  Vel X   Min: %f   Max: %f   [ kpc/h]\n", vx_min_g, vx_max_g);
  chprintf( "  Vel Y   Min: %f   Max: %f   [ kpc/h]\n", vy_min_g, vy_max_g);
  chprintf( "  Vel Z   Min: %f   Max: %f   [ kpc/h]\n", vz_min_g, vz_max_g);

  #endif


  free(dataset_buffer_px);
  free(dataset_buffer_py);
  free(dataset_buffer_pz);
  free(dataset_buffer_vx);
  free(dataset_buffer_vy);
  free(dataset_buffer_vz);

  #ifndef SINGLE_PARTICLE_MASS
  free(dataset_buffer_m);
  #endif
  #ifdef PARTICLE_IDS
  free(dataset_buffer_IDs);
  #endif
}
#endif






#endif
