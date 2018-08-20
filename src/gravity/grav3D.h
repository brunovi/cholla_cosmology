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

  // int n_ghost_potential;

  int n_cells;
  // int n_cells_potential;

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
  void Initialize( Real Lx, Real Ly, Real Lz, int nx_total, int ny_total, int nz_total, int nx_real, int ny_real, int nz_real, int dx_real, int dy_real, int dz_real);

  void Copy_global_parameters( Real dx_in, Real dy_in, Real dz_in, Real Lx, Real Ly, Real Lz );
  // /*! \fn void AllocateMemory(void)
  //  *  \brief Allocate memory for the d, m, E arrays. */
  // void AllocateMemory_CPU(void);
  // void AllocateMemory_GPU(void);
  //
  // // void CopyDensityFromHost( Grid3D G );
  // void Initialize_values_CPU();
  // void Initialize_values_GPU();
  //
  // #ifdef HDF5
  //   // /*! \fn void Write_GravHeader_HDF5(hid_t file_id)
  //   //  *  \brief Write the relevant header info to the HDF5 file. */
  //   void Write_GravHeader_HDF5(hid_t file_id);
  //
  //   /*! \fn void Write_GravFields_HDF5(hid_t file_id)
  //    *  \brief Write the grid to a file, at the current simulation time. */
  //   void Write_GravFields_HDF5(hid_t file_id);
  // #endif
  //
  // Real Get_Potential_PFFT( void );
  //
  // void Extrapolate_Potential( Real dt_0, Real dt_1, bool extrapolate );
  //
  // void Copy_Potential_To_Device( void );
  //
  // // void Copy_potential_to_Particles( void );
  //
  //
  // /*! \fn void Reset(void)
  //  *  \brief Reset the Griav3D class. */
  // void Reset(void);
  //
  // /*! \fn void FreeMemory(void)
  //  *  \brief Free the memory for the density array. */
  // void FreeMemory_CPU(void);
  // void FreeMemory_GPU(void);

};


#endif //GRAV3D_H
