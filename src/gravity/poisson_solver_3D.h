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


  void Initialize( Grav3D Grav );

};




#endif
#endif
