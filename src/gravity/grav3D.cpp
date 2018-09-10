#ifdef GRAVITY

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"../global.h"
#include"grav3D.h"
#include "../io.h"

// #include"self_gravity_functions.h"


Grav3D::Grav3D( void ){}


/*! \fn void Initialize(int nx_in, int ny_in, int nz_in)
 *  \brief Initialize the grid. */
void Grav3D::Initialize( Real Lx, Real Ly, Real Lz, int nx, int ny, int nz, int nx_real, int ny_real, int nz_real, Real dx_real, Real dy_real, Real dz_real, int n_ghost_pot_offset )
{

  Lbox_x = Lx;
  Lbox_y = Ly;
  Lbox_z = Lz;

  nx_total = nx;
  ny_total = ny;
  nz_total = nz;

  nx_local = nx_real;
  ny_local = ny_real;
  nz_local = nz_real;


  dx = dx_real;
  dy = dy_real;
  dz = dz_real;


  n_cells = nx_local*ny_local*nz_local;
  n_cells_potential = ( nx_local + 2*N_GHOST_POTENTIAL ) * ( ny_local + 2*N_GHOST_POTENTIAL ) * ( nz_local + 2*N_GHOST_POTENTIAL );

  INITIAL = true;
  dt_prev = 0;
  dt_now = 0;

  Gconst = 1;

  // Allocate memory
  AllocateMemory_CPU();

  Initialize_values_CPU();

  chprintf( "\nGravity Initialized: \n Lbox: %0.2f %0.2f %0.2f \n Local: %d %d %d \n Global: %d %d %d \n",
      Lbox_x, Lbox_y, Lbox_z, nx_local, ny_local, nz_local,   nx_total, ny_total, nz_total );


  chprintf( " dx:%f  dy:%f  dz:%f\n", dx, dy, dz );
  chprintf( " N ghost potential: %d\n", N_GHOST_POTENTIAL);
  chprintf( " N ghost offset: %d\n", n_ghost_pot_offset);
  // chprintf( "\n" );

}

void Grav3D::AllocateMemory_CPU(void)
{
  // allocate memory for the density and potential arrays
  F.density_h    = (Real *) malloc(n_cells*sizeof(Real));
  F.potential_h  = (Real *) malloc(n_cells_potential*sizeof(Real));
  F.potential_1_h  = (Real *) malloc(n_cells_potential*sizeof(Real));
}

void Grav3D::Initialize_values_CPU(void){

  for (int id=0; id<n_cells; id++){
    F.density_h[id] = 0;
  }

  for (int id_pot=0; id_pot<n_cells_potential; id_pot++){
    F.potential_h[id_pot] = 0;
    F.potential_1_h[id_pot] = 0;
  }
}

void Grav3D::FreeMemory_CPU(void)
{
  free(F.density_h);
  free(F.potential_h);
  free(F.potential_1_h);
}

#endif