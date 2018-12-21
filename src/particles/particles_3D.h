#ifdef PARTICLES

#ifndef PARTICLES_H
#define PARTICLES_H

#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <cstdlib>
#include <string.h>
#include"../global.h"
#include "../gravity/grav3D.h"

#ifdef GRAVITY_OMP
#include "../gravity/gravity_omp.h"
#endif

/*! \class Part3D
 *  \brief Class to create a set of particles in 3D space. */
class Particles_3D
{
  public:

  part_int_t n_local;

  part_int_t n_total;
  part_int_t n_total_0;

  int n_inside;
  int n_outside;

  Real dt;
  Real dt_max;
  Real t;

  Real C_cfl;

  #ifdef PARTICLE_IDS
  int_vector_t partIDs;
  #endif
  #ifndef SINGLE_PARTICLE_MASS
  real_vector_t mass;
  #endif
  real_vector_t pos_x;
  real_vector_t pos_y;
  real_vector_t pos_z;
  real_vector_t vel_x;
  real_vector_t vel_y;
  real_vector_t vel_z;
  real_vector_t grav_x;
  real_vector_t grav_y;
  real_vector_t grav_z;

  #ifdef REVERT_STEP
  #ifdef PARTICLE_IDS
  int_vector_t partIDs_0;
  #endif
  #ifndef SINGLE_PARTICLE_MASS
  real_vector_t mass_0;
  #endif
  real_vector_t pos_x_0;
  real_vector_t pos_y_0;
  real_vector_t pos_z_0;
  real_vector_t vel_x_0;
  real_vector_t vel_y_0;
  real_vector_t vel_z_0;
  real_vector_t grav_x_0;
  real_vector_t grav_y_0;
  real_vector_t grav_z_0;
  #endif

  bool INITIAL;

  #ifdef MPI_CHOLLA
  int_vector_t out_indxs_vec_x0;
  int_vector_t out_indxs_vec_x1;
  int_vector_t out_indxs_vec_y0;
  int_vector_t out_indxs_vec_y1;
  int_vector_t out_indxs_vec_z0;
  int_vector_t out_indxs_vec_z1;


  part_int_t n_transfer_x0;
  part_int_t n_transfer_x1;
  part_int_t n_transfer_y0;
  part_int_t n_transfer_y1;
  part_int_t n_transfer_z0;
  part_int_t n_transfer_z1;



  #endif
  bool TRANSFER_DENSITY_BOUNDARIES;
  bool TRANSFER_PARTICLES_BOUNDARIES;

  #ifdef COSMOLOGY
  Real current_z;
  Real current_a;
  #endif

  #ifdef SINGLE_PARTICLE_MASS
  Real particle_mass;
  #endif


  struct Grid
  {

    int nx_local, ny_local, nz_local;
    int nx_total, ny_total, nz_total;

    Real xMin, yMin, zMin;
    Real xMax, yMax, zMax;
    Real dx, dy, dz;

    Real domainMin_x, domainMax_x;
    Real domainMin_y, domainMax_y;
    Real domainMin_z, domainMax_z;

    int n_ghost_particles_grid;
    int n_cells;

    Real *density;
    Real *gravity_x;
    Real *gravity_y;
    Real *gravity_z;


  } G;

  Particles_3D(void);

  void Initialize( struct parameters P, Grav3D &Grav, Real xblocal, Real yblocal, Real zblocal, Real xbound, Real ybound, Real zbound, Real xdglobal, Real ydglobal, Real zdglobal );

  void AllocateMemory_CPU( void );
  void FreeMemory_CPU( void );

  void Initialize_values_CPU( void );
  void Initialize_Sphere( void );
  void Initialize_Uniform_Grid( void );

  #ifdef MPI_CHOLLA
  void Clear_Vectors_For_Transfers( void );
  // void Select_Particles_to_Transfer( void );
  void Select_Particles_to_Transfer( int dir );
  void Load_Particles_to_Buffer( int direction, int side, int buffer_start, Real *send_buffer, int MAX_PARTICLES_IN_BUFFER , bool secondary );
  void Load_Particles_to_Buffer_X( void );
  void Load_Particles_to_Buffer_Y( void );
  void Load_Particles_to_Buffer_Z( void );
  void Unload_Particles_from_Buffer( int direction, int side, int buffer_start, Real *recv_buffer, Real *send_buffer_y0, Real *send_buffer_y1, Real *send_buffer_z0, Real *send_buffer_z1, int buffer_start_y, int buffer_start_z, Real *send_buffer_y0_second, Real *send_buffer_y1_second, Real *send_buffer_z0_second, Real *send_buffer_z1_second, int max_particles, bool secondary  );
  void Add_Particle_To_Buffer( Real *buffer, int buffer_start, int max_particles, Real *buffer_secondary, int buffer_start_secondary, Real pId, Real pMass, Real pPos_x, Real pPos_y, Real pPos_z, Real pVel_x, Real pVel_y, Real pVel_z, bool secondary );
  void Add_Particle_To_Vectors( Real pId, Real pMass, Real pPos_x, Real pPos_y, Real pPos_z, Real pVel_x, Real pVel_y, Real pVel_z );
  void Remove_Transfered_Particles( void );
  void Clear_Particles_For_Transfer( void );
  #endif

  void Reset( void );


};






#endif
#endif //PARTICLES_H
