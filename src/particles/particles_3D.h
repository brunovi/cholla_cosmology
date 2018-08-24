#ifdef PARTICLES

#ifndef PARTICLES_H
#define PARTICLES_H

#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <cstdlib>
#include"../global.h"
#include "../gravity/grav3D.h"




#include <vector>


typedef std::vector<Real> real_vector_t;
typedef std::vector<part_int_t> int_vector_t;


/*! \class Part3D
 *  \brief Class to create a set of particles in 3D space. */
class Particles_3D
{
  public:

  part_int_t n_local;

  part_int_t n_total;

  int n_inside;
  int n_outside;

  Real dt;
  Real dt_max;
  Real t;

  int_vector_t partIDs;
  real_vector_t mass;
  real_vector_t pos_x;
  real_vector_t pos_y;
  real_vector_t pos_z;
  real_vector_t vel_x;
  real_vector_t vel_y;
  real_vector_t vel_z;
  real_vector_t grav_x;
  real_vector_t grav_y;
  real_vector_t grav_z;


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

    /*! \var density
     *  \brief Array containing the density of each cell in the grid */
    Real *density_h;
    Real *gravity_x;
    Real *gravity_y;
    Real *gravity_z;


  } G;

  Particles_3D(void);

  void Initialize( Grav3D &Grav, Real xblocal, Real yblocal, Real zblocal, Real xbound, Real ybound, Real zbound, Real xdglobal, Real ydglobal, Real zdglobal );



};






#endif
#endif //PARTICLES_H
