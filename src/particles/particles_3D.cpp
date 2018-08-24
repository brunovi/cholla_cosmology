#ifdef PARTICLES


#include "particles_3D.h"


Particles_3D::Particles_3D( void ){}


void Particles_3D::Initialize( Grid3D &Grid_hydro ){

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

  G.nx_local = Grid_hydro.Grav.nx_local;
  G.ny_local = Grid_hydro.Grav.ny_local;
  G.nz_local = Grid_hydro.Grav.nz_local;
  G.nx_total = Grid_hydro.Grav.nx_total;
  G.ny_total = Grid_hydro.Grav.ny_total;
  G.nz_total = Grid_hydro.Grav.nz_total;

  G.dx = Grid_hydro.H.dx;
  G.dy = Grid_hydro.H.dy;
  G.dz = Grid_hydro.H.dz;

  G.xMin = Grid_hydro.H.xblocal;
  G.yMin = Grid_hydro.H.yblocal;
  G.zMin = Grid_hydro.H.zblocal;

  G.xMax = G.xMin + G.nx_local*G.dx;
  G.yMax = G.yMin + G.ny_local*G.dy;
  G.zMax = G.zMin + G.nz_local*G.dz;

  G.domainMin_x = Grid_hydro.H.xbound;
  G.domainMax_x = Grid_hydro.H.xbound + Grid_hydro.H.xdglobal;
  G.domainMin_y = Grid_hydro.H.ybound;
  G.domainMax_y = Grid_hydro.H.ybound + Grid_hydro.H.ydglobal;
  G.domainMin_z = Grid_hydro.H.zbound;
  G.domainMax_z = Grid_hydro.H.zbound + Grid_hydro.H.zdglobal;

  chprintf("Particles Initialized: \n n_local: %lu \n", n_local );
  chprintf(" xDomain_local:  [%.4f %.4f ] [%.4f %.4f ] [%.4f %.4f ]\n", G.xMin, G.xMax, G.yMin, G.yMax, G.zMin, G.zMax );
  chprintf(" xDomain_global: [%.4f %.4f ] [%.4f %.4f ] [%.4f %.4f ]\n  ", G.domainMin_x, G.domainMax_x, G.domainMin_y, G.domainMax_y, G.domainMin_z, G.domainMax_z);

}


#endif
