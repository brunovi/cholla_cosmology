#ifndef GRAV3D_H
#define GRAV3D_H

#include<stdio.h>
#include"../global.h"


#ifdef HDF5
#include<hdf5.h>
#endif



/*! \class Grid3D
 *  \brief Class to create a 3D grid of cells. */
class Grav3D
{
  public:

  Real Lbox_x;
  Real Lbox_y;
  Real Lbox_z;
  /*! \var nx
  *  \brief Total number of cells in the x-dimension */
  int nx_total;
  /*! \var ny
  *  \brief Total number of cells in the y-dimension */
  int ny_total;
  /*! \var nz
  *  \brief Total number of cells in the z-dimension */
  int nz_total;

  /*! \var nx_local
  *  \brief Local number of cells in the x-dimension */
  int nx_local;
  /*! \var ny_local
  *  \brief Local number of cells in the y-dimension */
  int ny_local;
  /*! \var nz_local
  *  \brief Local number of cells in the z-dimension */
  int nz_local;

  /*! \var dx
  *  \brief x-width of cells */
  Real dx;
  /*! \var dy
  *  \brief y-width of cells */
  Real dy;
  /*! \var dz
  *  \brief z-width of cells */
  Real dz;

  #ifdef COSMOLOGY
  Real current_a;
  #endif

  Real dens_avrg ;


  int n_cells;
  int n_cells_potential;


  bool INITIAL;

  Real dt_prev;
  Real dt_now;

  Real Gconst;

  struct Fields
  {
    /*! \var density
     *  \brief Array containing the density of each cell in the grid */
    Real *density_h;

    /*! \var potential_h
     *  \brief Array containing the gravitational potential of each cell in the grid */
    Real *potential_h;

    /*! \var potential_h
     *  \brief Array containing the gravitational potential of each cell in the grid */
    Real *potential_1_h;


  } F;

  /*! \fn Grav3D(void)
   *  \brief Constructor for the gravity class */
  Grav3D(void);

  /*! \fn void Initialize(int nx_in, int ny_in, int nz_in)
   *  \brief Initialize the grid. */
  void Initialize( Real Lx, Real Ly, Real Lz, int nx_total, int ny_total, int nz_total, int nx_real, int ny_real, int nz_real, Real dx_real, Real dy_real, Real dz_real, int n_ghost_pot_offset);

  // /*! \fn void AllocateMemory(void)
  void AllocateMemory_CPU(void);

  void Initialize_values_CPU();
  void FreeMemory_CPU(void);
};


#endif //GRAV3D_H
