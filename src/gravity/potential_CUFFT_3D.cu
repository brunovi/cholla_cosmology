#ifdef GRAVITY
#ifdef POTENTIAL_CUFFT

#include "potential_CUFFT_3D.h"


Potential_CUFFT_3D::Potential_CUFFT_3D( void ){}

void Potential_CUFFT_3D::Initialize( Grav3D Grav){

  Lbox_x = Grav.Lbox_x;
  Lbox_y = Grav.Lbox_y;
  Lbox_z = Grav.Lbox_z;

  nx_total = Grav.nx_total;
  ny_total = Grav.ny_total;
  nz_total = Grav.nz_total;

  nx_local = Grav.nx_local;
  ny_local = Grav.ny_local;
  nz_local = Grav.nz_local;

  dx = Grav.dx;
  dy = Grav.dy;
  dz = Grav.dz;

  chprintf( " Using Poisson Solver: CUFFT\n");
chprintf( "  CUFFT: L[ %f %f %f ] N[ %d %d %d ] dx[ %f %f %f ]\n", Lbox_x, Lbox_y, Lbox_z, nx_local, ny_local, nz_local, dx, dy, dz );



}




#endif //POTENTIAL_CUFFT
#endif //GRAVITY
