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
void Grav3D::Initialize( Real x_min, Real y_min, Real z_min, Real Lx, Real Ly, Real Lz, int nx, int ny, int nz, int nx_real, int ny_real, int nz_real, Real dx_real, Real dy_real, Real dz_real, int n_ghost_pot_offset )
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

  xMin = x_min;
  yMin = y_min;
  zMin = z_min;


  n_cells = nx_local*ny_local*nz_local;
  n_cells_potential = ( nx_local + 2*N_GHOST_POTENTIAL ) * ( ny_local + 2*N_GHOST_POTENTIAL ) * ( nz_local + 2*N_GHOST_POTENTIAL );

  INITIAL = true;
  dt_prev = 0;
  dt_now = 0;

  #ifdef COSMOLOGY
  current_a = 0;
  #endif

  dens_avrg = 0;

  Gconst = 1;

  TRANSFER_POTENTIAL_BOUNDARIES = false;

  // Allocate memory
  AllocateMemory_CPU();

  Initialize_values_CPU();

  chprintf( "\nGravity Initialized: \n Lbox: %0.2f %0.2f %0.2f \n Local: %d %d %d \n Global: %d %d %d \n",
      Lbox_x, Lbox_y, Lbox_z, nx_local, ny_local, nz_local,   nx_total, ny_total, nz_total );


  chprintf( " dx:%f  dy:%f  dz:%f\n", dx, dy, dz );
  chprintf( " N ghost potential: %d\n", N_GHOST_POTENTIAL);
  chprintf( " N ghost offset: %d\n", n_ghost_pot_offset);
  // std::cout << xMin << std::endl;
  // chprintf( "\n" );

}

void Grav3D::AllocateMemory_CPU(void)
{
  // allocate memory for the density and potential arrays
  F.density_h    = (Real *) malloc(n_cells*sizeof(Real));
  F.potential_h  = (Real *) malloc(n_cells_potential*sizeof(Real));
  F.potential_1_h  = (Real *) malloc(n_cells_potential*sizeof(Real));

  #ifdef GRAVITY_CORRECTOR
  F.gravity_x_h  = (Real *) malloc(n_cells*sizeof(Real));
  F.gravity_y_h  = (Real *) malloc(n_cells*sizeof(Real));
  F.gravity_z_h  = (Real *) malloc(n_cells*sizeof(Real));
  F.gravity_x_h_prev  = (Real *) malloc(n_cells*sizeof(Real));
  F.gravity_y_h_prev  = (Real *) malloc(n_cells*sizeof(Real));
  F.gravity_z_h_prev  = (Real *) malloc(n_cells*sizeof(Real));
  #endif

  #ifdef EXTRA_FIELD
  F.extra_field  = (Real *) malloc(n_cells*sizeof(Real));
  #endif

}

void Grav3D::Initialize_values_CPU(void){

  for (int id=0; id<n_cells; id++){
    F.density_h[id] = 0;
    #ifdef GRAVITY_CORRECTOR
    F.gravity_x_h[id] = 0;
    F.gravity_y_h[id] = 0;
    F.gravity_z_h[id] = 0;
    F.gravity_x_h_prev[id] = 0;
    F.gravity_y_h_prev[id] = 0;
    F.gravity_z_h_prev[id] = 0;
    #endif
    #ifdef EXTRA_FIELD
    F.extra_field[id] = 0;
    #endif
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

  #ifdef GRAVITY_CORRECTOR
  free(F.gravity_x_h);
  free(F.gravity_y_h);
  free(F.gravity_z_h);
  free(F.gravity_x_h_prev);
  free(F.gravity_y_h_prev);
  free(F.gravity_z_h_prev);
  #endif

  #ifdef EXTRA_FIELD
  free(F.extra_field);
  #endif
}

#endif
