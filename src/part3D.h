//
// #ifndef PARTICLES_H
// #define PARTICLES_H
//
// #include<stdio.h>
// #include"global.h"
// #include "grid3D.h"
// #include "io.h"
//
// #ifdef HDF5
// #include<hdf5.h>
// #endif
//
// #include <vector>
//
//
// typedef std::vector<Real> ch_vector_t;
// typedef std::vector<partID_t> ids_vector_t;
//
//
// /*! \class Grid3D
//  *  \brief Class to create a 3D grid of cells. */
// class Part3D
// {
//   public:
//
//   /*! \var nx
//   *  \brief Total number of cells in the x-dimension */
//   partID_t n_local;
//   /*! \var ny
//   *  \brief Total number of cells in the y-dimension */
//   partID_t n_total;
//   /*! \var nz
//   *  \brief Total number of cells in the z-dimension */
//   int n_ghost;
//   int n_inside;
//   int n_outside;
//
//   Real dt;
//   Real dt_max;
//   Real t;
//
//   ids_vector_t partIDs;
//   ch_vector_t mass;
//   ch_vector_t pos_x;
//   ch_vector_t pos_y;
//   ch_vector_t pos_z;
//   ch_vector_t vel_x;
//   ch_vector_t vel_y;
//   ch_vector_t vel_z;
//   ch_vector_t grav_x;
//   ch_vector_t grav_y;
//   ch_vector_t grav_z;
//
//   ids_vector_t out_indxs_vec_x0;
//   ids_vector_t out_indxs_vec_x1;
//   ids_vector_t out_indxs_vec_y0;
//   ids_vector_t out_indxs_vec_y1;
//   ids_vector_t out_indxs_vec_z0;
//   ids_vector_t out_indxs_vec_z1;
//
//   struct Grid
//   {
//
//     int nx, ny, nz;
//     int nx_total, ny_total, nz_total;
//     int grid_ghost;
//     int grid_ghost_pot;
//     Real xMin, yMin, zMin;
//     Real xMax, yMax, zMax;
//     Real dx, dy, dz;
//
//     Real domainMin_x, domainMax_x;
//     Real domainMin_y, domainMax_y;
//     Real domainMin_z, domainMax_z;
//
//     /*! \var density
//      *  \brief Array containing the density of each cell in the grid */
//     Real *density_h;
//     Real *gravity_x;
//     Real *gravity_y;
//     Real *gravity_z;
//
//     Real *potential;
//     Real *potential_cosmo;
//
//     // /*! \var potential_h
//     //  *  \brief Array containing the gravitational potential of each cell in the grid */
//     // Real *potential_h;
//     //
//     // /*! \var potential_h
//     //  *  \brief Array containing the gravitational potential of each cell in the grid */
//     // Real *potential_1_h;
//     //
//     // /*! \var potential_h
//     //  *  \brief Array containing the gravitational potential in the GPU of each cell in the grid */
//     // Real *potential_d;
//
//
//   } G;
//
//   // struct Prop
//   // {
//   //   /*! \var density
//   //    *  \brief Array containing the density of each cell in the grid */
//   //   ch_vector_t pos_x;
//   //   ch_vector_t pos_y;
//   //   ch_vector_t pos_z;
//   //   ch_vector_t vel_x;
//   //   ch_vector_t vel_y;
//   //   ch_vector_t vel_z;
//   //
//   // } Prop;
//
//   /*! \fn Grav3D(void)
//    *  \brief Constructor for the gravity class */
//   Part3D(void);
//
//   /*! \fn void Initialize(int nx_in, int ny_in, int nz_in)
//    *  \brief Initialize the grid. */
//   void Initialize( partID_t n_local, partID_t n_total, Grid3D &Grid_hydro );
//
//   /*! \fn void AllocateMemory(void)
//    *  \brief Allocate memory for the d, m, E arrays. */
//   void AllocateMemory_CPU(void);
//   void Allocate_Buffers_For_Transfers( void );
//
//   // void CopyDensityFromHost( Grid3D G );
//   void Initialize_values_CPU( Grid3D &Grid_hydro );
//
//   void Load_PartData(struct parameters P);
//
//
//   #ifdef HDF5
//     // /*! \fn void Write_GravHeader_HDF5(hid_t file_id)
//     //  *  \brief Write the relevant header info to the HDF5 file. */
//     void Write_PartHeader_HDF5(hid_t file_id);
//
//     /*! \fn void Write_GravFields_HDF5(hid_t file_id)
//      *  \brief Write the grid to a file, at the current simulation time. */
//     void Write_PartData_HDF5(hid_t file_id);
//     void Load_PartData_HDF5(hid_t file_id, int nfile, bool print_load);
//   #endif
//
//
//   void Get_Indexes_CIC( Real pos_x, Real pos_y, Real pos_z, int &indx_x, int &indx_y, int &indx_z );
//   void Clear_Density( void );
//   void Get_Density_CIC( partID_t p_indx_start, partID_t p_indx_end );
//   void Copy_potential_from_Grid( Grid3D &Grid );
//   void Copy_potential_to_Grid( Grid3D &Grid );
//   void Get_Potential_Gradient( bool cosmological, int grid_indx_start, int grid_indx_end );
//   #ifdef COSMOLOGY
//   // void Get_Potential_Gradient_cosmo( void );
//   void potential_proper_to_cosmological( void );
//   void Advance_Particle_Step_LeapFrog_cosmo( partID_t pIndx, bool only_vel );
//   #endif
//   void Get_Gravity_CIC( partID_t p_indx_start, partID_t p_indx_end );
//   void Add_Density_to_Grid( Grid3D &Grid, bool clearPrevious );
//
//   void Advance_Particle_Step_LeapFrog( partID_t pIndx, Real dt, bool only_vel );
//   void Advance_Step_LeapFrog( Real dt, bool only_vel, partID_t p_indx_start, partID_t p_indx_end );
//
//   void Select_Particles_to_Transfer( void );
//   void Load_Particles_to_Buffer( int direction, int side, ids_vector_t &out_indxs_vec );
//   void Load_Particles_to_Buffer_X( void );
//   void Load_Particles_to_Buffer_Y( void );
//   void Load_Particles_to_Buffer_Z( void );
//   void Finish_Unload_from_Transfers( void );
//
//   void Unload_Particles_from_Buffer( int direction, int side );
//
//   void Add_Particle_To_Vectors( Real pId, Real pMass, Real pPos_x, Real pPos_y, Real pPos_z,
//                               Real pVel_x, Real pVel_y, Real pVel_z );
//
//   void Remove_Particles( void );
//
//
//
//
//   // Real Get_Potential_PFFT( void );
//
//   // void Extrapolate_Potential( Real dt_0, Real dt_1 );
//
//   // void Copy_Potential_To_Device( void );
//
//
//   /*! \fn void Reset(void)
//    *  \brief Reset the Griav3D class. */
//   void Reset(void);
//
//   /*! \fn void FreeMemory(void)
//    *  \brief Free the memory for the density array. */
//   void FreeMemory_CPU(void);
//   // void FreeMemory_GPU(void);
//
// };
//
//
// #endif //PARTICLES_H
