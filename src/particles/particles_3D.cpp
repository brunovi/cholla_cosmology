#ifdef PARTICLES


#include "particles_3D.h"
#include "../io.h"


Particles_3D::Particles_3D( void ){}


void Particles_3D::Initialize( Grav3D &Grav, Real xblocal, Real yblocal, Real zblocal, Real xbound, Real ybound, Real zbound, Real xdglobal, Real ydglobal, Real zdglobal ){

  chprintf( "\nInitializing Particles...\n");
  n_local = 0;
  n_total = 0;

  dt = 0.0;
  t = 0.0;
  dt_max = 0.01;

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

  G.nx_local = Grav.nx_local;
  G.ny_local = Grav.ny_local;
  G.nz_local = Grav.nz_local;
  G.nx_total = Grav.nx_total;
  G.ny_total = Grav.ny_total;
  G.nz_total = Grav.nz_total;

  G.dx = Grav.dx;
  G.dy = Grav.dy;
  G.dz = Grav.dz;

  G.xMin = xblocal;
  G.yMin = yblocal;
  G.zMin = zblocal;

  G.xMax = G.xMin + G.nx_local*G.dx;
  G.yMax = G.yMin + G.ny_local*G.dy;
  G.zMax = G.zMin + G.nz_local*G.dz;

  G.domainMin_x = xbound;
  G.domainMax_x = xbound + xdglobal;
  G.domainMin_y = ybound;
  G.domainMax_y = ybound + ydglobal;
  G.domainMin_z = zbound;
  G.domainMax_z = zbound + zdglobal;

  chprintf("Particles Initialized: \n n_local: %lu \n", n_local );
  chprintf(" xDomain_local:  [%.4f %.4f ] [%.4f %.4f ] [%.4f %.4f ]\n", G.xMin, G.xMax, G.yMin, G.yMax, G.zMin, G.zMax );
  chprintf(" xDomain_global: [%.4f %.4f ] [%.4f %.4f ] [%.4f %.4f ]\n", G.domainMin_x, G.domainMax_x, G.domainMin_y, G.domainMax_y, G.domainMin_z, G.domainMax_z);
  chprintf(" dx: %f  %f  %f\n", G.dx, G.dy, G.dz );
}


#endif
