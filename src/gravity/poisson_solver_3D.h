#ifdef GRAVITY

#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "../global.h"
#include "grav3D.h"


class Poisson_Solver_3D{

  public:

  Real Lbox_x;
  Real Lbox_y;
  Real Lbox_z;

  int nx_total;
  int ny_total;
  int nz_total;

  int nx_local;
  int ny_local;
  int nz_local;

  Real dx;
  Real dy;
  Real dz;
  int n_cells_total;
  int n_cells_local;

  // Poisson_Solver_3D( void );
  void Initialize( Grav3D Grav );

  void Get_Potential( Grav3D Grav );



};






#endif
#endif
