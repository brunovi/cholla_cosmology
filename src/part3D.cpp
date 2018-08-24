// #ifdef PARTICLES
//
// #include<stdio.h>
// #include<stdlib.h>
// #include"global.h"
// #include <math.h>
// #include <cstdlib>
// #include "part3D.h"
// #include "random_routines.h"
//
// #ifdef COSMOLOGY
// #include"cosmology.h"
// #endif
//
// Real Get_Random_Uniform( Real a, Real b){
//   Real val_int = (Real) rand();
//   return val_int/RAND_MAX * ( b - a ) + a;
// }
//
//
// Part3D::Part3D( void ){}
//
//
// void Part3D::Initialize( partID_t n_local_0, partID_t n_total_0, Grid3D &Grid_hydro ){
//
//   n_local = n_local_0;
//   n_total = n_total_0;
//
//   dt = 0.0;
//   t = 0.0;
//   dt_max = 0.01;
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
//
//
//   out_indxs_vec_x0.clear();
//   out_indxs_vec_x1.clear();
//   out_indxs_vec_y0.clear();
//   out_indxs_vec_y1.clear();
//   out_indxs_vec_z0.clear();
//   out_indxs_vec_z1.clear();
//
//
//   G.nx = Grid_hydro.Grav.nx_local;
//   G.ny = Grid_hydro.Grav.ny_local;
//   G.nz = Grid_hydro.Grav.nz_local;
//   G.nx_total = Grid_hydro.Grav.nx_total;
//   G.ny_total = Grid_hydro.Grav.ny_total;
//   G.nz_total = Grid_hydro.Grav.nz_total;
//
//   G.grid_ghost = 1;
//   G.grid_ghost_pot = 2;
//
//   G.dx = Grid_hydro.H.dx;
//   G.dy = Grid_hydro.H.dy;
//   G.dz = Grid_hydro.H.dz;
//
//   G.xMin = Grid_hydro.H.xblocal;
//   G.yMin = Grid_hydro.H.yblocal;
//   G.zMin = Grid_hydro.H.zblocal;
//
//   G.xMax = G.xMin + G.nx*G.dx;
//   G.yMax = G.yMin + G.ny*G.dy;
//   G.zMax = G.zMin + G.nz*G.dz;
//
//   G.domainMin_x = Grid_hydro.H.xbound;
//   G.domainMax_x = Grid_hydro.H.xbound + Grid_hydro.H.xdglobal;
//   G.domainMin_y = Grid_hydro.H.ybound;
//   G.domainMax_y = Grid_hydro.H.ybound + Grid_hydro.H.ydglobal;
//   G.domainMin_z = Grid_hydro.H.zbound;
//   G.domainMax_z = Grid_hydro.H.zbound + Grid_hydro.H.zdglobal;
//
//   #ifdef COSMO_UNITS
//   // G.dx /= r_0;
//   // G.dy /= r_0;
//   // G.dz /= r_0;
//   // G.xMin /= r_0;
//   // G.yMin /= r_0;
//   // G.zMin /= r_0;
//   // G.xMax /= r_0;
//   // G.yMax /= r_0;
//   // G.zMax /= r_0;
//   // G.domainMin_x /= r_0;
//   // G.domainMax_x /= r_0;
//   // G.domainMin_y /= r_0;
//   // G.domainMax_y /= r_0;
//   // G.domainMin_z /= r_0;
//   // G.domainMax_z /= r_0;
//   #endif
//
//   AllocateMemory_CPU();
//   Allocate_Buffers_For_Transfers();
//   Initialize_values_CPU( Grid_hydro );
//
//   Grid_hydro.Grav.F.density_particles = G.density_h;
//   Grid_hydro.Grav.F.gravX_particles = G.gravity_x;
//   Grid_hydro.Grav.F.gravY_particles = G.gravity_y;
//   Grid_hydro.Grav.F.gravZ_particles = G.gravity_z;
//   Grid_hydro.Grav.F.potential_particles = G.potential;
//
//   chprintf("Particles Initialized: \n n_local: %lu \n", n_local );
//   chprintf(" xDomain_local:  [%.4f %.4f ]  dx: %.4f\n", G.xMin, G.xMax, G.dx );
//   chprintf(" xDomain_global: [%.4f %.4f ]  dx: %.4f\n", G.domainMin_x, G.domainMax_x);
//
// }
//
// void Part3D::AllocateMemory_CPU( void ){
//   int nGHST = G.grid_ghost;
//   int nGHST_pot = G.grid_ghost_pot;
//   int n_cells = (G.nx+2*nGHST) * (G.ny+2*nGHST) * (G.nz+2*nGHST);
//   int n_cells_p = (G.nx+2*nGHST_pot) * (G.ny+2*nGHST_pot) * (G.nz+2*nGHST_pot);
//   G.density_h = (Real *) malloc(n_cells*sizeof(Real));
//   G.gravity_x = (Real *) malloc(n_cells*sizeof(Real));
//   G.gravity_y = (Real *) malloc(n_cells*sizeof(Real));
//   G.gravity_z = (Real *) malloc(n_cells*sizeof(Real));
//   G.potential = (Real *) malloc(n_cells_p*sizeof(Real));
//   G.potential_cosmo = (Real *) malloc(n_cells_p*sizeof(Real));
//
// }
//
// void Part3D::Initialize_values_CPU( Grid3D &Grid_hydro )
// {
//   int i, j, k, id_dens, id_pot;
//   partID_t pId;
//   Real x_pos, y_pos, z_pos;
//   int nGHST = G.grid_ghost;
//   int nGHST_pot = G.grid_ghost_pot;
//   int n_cells = (G.nx+2*nGHST) * (G.ny+2*nGHST) * (G.nz+2*nGHST);
//   int n_cells_pot = (G.nx+2*nGHST_pot) * (G.ny+2*nGHST_pot) * (G.nz+2*nGHST_pot);
//
//   for( i=0; i<n_cells; i++ ){
//     G.density_h[i] = 0;
//     G.gravity_x[i] = 0;
//     G.gravity_y[i] = 0;
//     G.gravity_z[i] = 0;
//   }
//
//   for( i=0; i<n_cells_pot; i++ ){
//     G.potential[i] = 0;
//   }
//
//
//   Real center_x, center_y, center_z, radius, sphereR;
//   center_x = 500;
//   center_y = 500;
//   center_z = 500;
//   // sphereR = 115000.0/5;
//   sphereR = 200;
//   //
//   //
//   // n_local = G.nx*G.ny*G.nz;
//   // Real pPos_x, pPos_y, pPos_z;
//   // chprintf (" nParticles_local: %ld\n", n_local);
//   //
//   // // Real rho_start = 1;
//   // // Real M_sphere = 4./3 * M_PI* rho_start * sphereR*sphereR*sphereR;
//   // // Real Mparticle = M_sphere / (8*n_local);
//   //
//   // Real r_halo = 200;
//   // Real M_halo = 1e12 / 2;
//   // Real rho_halo = M_halo / (4./3 * M_PI * r_halo*r_halo*r_halo);
//   // Real Mparticle = M_halo / (8*n_local);
//
//   // partID_t counter = 0;
//   // while (counter < n_local) {
//   //   pPos_x = rand_real( G.xMin, G.xMax );
//   //   pPos_y = rand_real( G.yMin, G.yMax );
//   //   pPos_z = rand_real( G.zMin, G.zMax );
//   //   radius = sqrt( (pPos_x-center_x)*(pPos_x-center_x) + (pPos_y-center_y)*(pPos_y-center_y) + (pPos_z-center_z)*(pPos_z-center_z) );
//   //   if ( radius > sphereR ) continue;
//   //   partIDs.push_back( (partID_t) procID*n_local + counter );
//   //   mass.push_back( Mparticle );
//   //   pos_x.push_back( pPos_x );
//   //   pos_y.push_back( pPos_y );
//   //   pos_z.push_back( pPos_z);
//   //   vel_x.push_back( 0.0 );
//   //   vel_y.push_back( 0.0 );
//   //   vel_z.push_back( 0.0 );
//   //   grav_x.push_back( 0.0 );
//   //   grav_y.push_back( 0.0 );
//   //   grav_z.push_back( 0.0 );
//   //   counter += 1;
//   // }
//
//
//   // Real pPos_x_loc, pPos_y_loc, pPos_z_loc;
//   // for (k=0; k<G.nz; k++) {
//   //   for (j=0; j<G.ny; j++) {
//   //     for (i=0; i<G.nx; i++) {
//   //       id_pot = (i+nGHST_pot) + (j+nGHST_pot)*(G.nx+2*nGHST_pot) + (k+nGHST_pot)*(G.nx+2*nGHST_pot)*(G.ny+2*nGHST_pot);
//   //       pId = (i) + (j)*G.nx + (k)*G.nx*G.ny;
//   //       pPos_x = G.xMin + i*G.dx + 0.5*G.dx;
//   //       pPos_y = G.yMin + j*G.dy + 0.5*G.dy;
//   //       pPos_z = G.zMin + k*G.dz + 0.5*G.dz;
//         // if ( k != 0 ) continue;
//         // if ( j != G.ny/2 ) continue;
//         // if ( i != G.nx/2 ) continue;
//         // partIDs.push_back( (partID_t) procID );
//         // mass.push_back( 1.0 );
//         // pos_x.push_back( pPos_x );
//         // pos_y.push_back( pPos_y );
//         // pos_z.push_back( pPos_z);
//         // vel_x.push_back( 0.0 );
//         // vel_y.push_back( 0.0 );
//         // vel_z.push_back( 0.0 );
//         // grav_x.push_back( 0.0 );
//         // grav_y.push_back( 0.0 );
//         // grav_z.push_back( 0.0 );
//         // n_local += 1;
//         // radius = sqrt( (pPos_x-center_x)*(pPos_x-center_x) + (pPos_y-center_y)*(pPos_y-center_y) + (pPos_z-center_z)*(pPos_z-center_z) );
//         // if ( radius < sphereR ){
//           // pPos_x += Get_Random_Uniform( -0.2*G.dx, 0.2*G.dx );
//           // pPos_y += Get_Random_Uniform( -0.2*G.dy, 0.2*G.dy );
//           // pPos_z += Get_Random_Uniform( -0.2*G.dz, 0.2*G.dz );
//           // pPos_x += 0.2*G.dx;
//           // pPos_y += 0.2*G.dy;
//           // pPos_z += 0.2*G.dz;
//           // partIDs.push_back(pId + procID*n_local);
//           // mass.push_back( 1.0 );
//           // pos_x.push_back( pPos_x );
//           // pos_y.push_back( pPos_y );
//           // pos_z.push_back( pPos_z);
//           // vel_x.push_back( 0.0 );
//           // vel_y.push_back( 0.0 );
//           // vel_z.push_back( 0.0 );
//           // grav_x.push_back( 0.0 );
//           // grav_y.push_back( 0.0 );
//           // grav_z.push_back( 0.0 );
//           // n_local += 1;
//         // }
//         // if ( i==0 )    vel_x[pId] = -1;
//         // if ( i==G.nx-1 ) vel_x[pId] =  1;
//         // if ( j==0 )    vel_y[pId] = -1;
//         // if ( j==G.ny-1 ) vel_y[pId] =  1;
//         // if ( k==0 )    vel_z[pId] = -1;
//         // if ( k==G.nz-1 ) vel_z[pId] =  1;
//         // G.potential[id_pot] = 1;
//   //     }
//   //   }
//   // }
//
// }
//
//
//
// void Part3D::FreeMemory_CPU(void)
// {
//  free(G.density_h);
//  free(G.gravity_x);
//  free(G.gravity_y);
//  free(G.gravity_z);
//  free(G.potential);
//  free(G.potential_cosmo);
//  // free(F.potential_1_h);
//
//  free(send_buffer_x0_particles);
//  free(send_buffer_x1_particles);
//  free(recv_buffer_x0_particles);
//  free(recv_buffer_x1_particles);
//  free(send_buffer_y0_particles);
//  free(send_buffer_y1_particles);
//  free(recv_buffer_y0_particles);
//  free(recv_buffer_y1_particles);
//  free(send_buffer_z0_particles);
//  free(send_buffer_z1_particles);
//  free(recv_buffer_z0_particles);
//  free(recv_buffer_z1_particles);
// }
//
// // #ifdef HDF5
// // #include<hdf5.h>
// // #endif
// //
// // #include <vector>
// //
// // std::vector<Real> ch_vector;
// //
// //
// // /*! \class Grid3D
// //  *  \brief Class to create a 3D grid of cells. */
// // class Particles
// // {
// //   public:
// //
// //   /*! \var nx
// //   *  \brief Total number of cells in the x-dimension */
// //   int n_local;
// //   /*! \var ny
// //   *  \brief Total number of cells in the y-dimension */
// //   int n_total;
// //   /*! \var nz
// //   *  \brief Total number of cells in the z-dimension */
// //   int n_ghost;
// //   int n_inside;
// //   int n_outside;
// //
// //
// //   struct Prop
// //   {
// //     /*! \var density
// //      *  \brief Array containing the density of each cell in the grid */
// //     ch_vector *pos_x;
// //     ch_vector *pos_y;
// //     ch_vector *pos_z;
// //     ch_vector *vel_x;
// //     ch_vector *vel_y;
// //     ch_vector *vel_z;
// //
// //   } Prop;
// //
// //   /*! \fn Grav3D(void)
// //    *  \brief Constructor for the gravity class */
// //   Particles(void);
// //
// //   /*! \fn void Initialize(int nx_in, int ny_in, int nz_in)
// //    *  \brief Initialize the grid. */
// //   void Initialize( int n_local, int n_total );
// //
// //   /*! \fn void AllocateMemory(void)
// //    *  \brief Allocate memory for the d, m, E arrays. */
// //   // void AllocateMemory_CPU(void);
// //
// //
// //   // void CopyDensityFromHost( Grid3D G );
// //   void Initialize_values_CPU();
// //
// //
// //   #ifdef HDF5
// //     // /*! \fn void Write_GravHeader_HDF5(hid_t file_id)
// //     //  *  \brief Write the relevant header info to the HDF5 file. */
// //     void Write_PartHeader_HDF5(hid_t file_id);
// //
// //     /*! \fn void Write_GravFields_HDF5(hid_t file_id)
// //      *  \brief Write the grid to a file, at the current simulation time. */
// //     void Write_PartFields_HDF5(hid_t file_id);
// //   #endif
// //
// //   // Real Get_Potential_PFFT( void );
// //
// //   // void Extrapolate_Potential( Real dt_0, Real dt_1 );
// //
// //   // void Copy_Potential_To_Device( void );
// //
// //
// //   /*! \fn void Reset(void)
// //    *  \brief Reset the Griav3D class. */
// //   void Reset(void);
// //
// //   /*! \fn void FreeMemory(void)
// //    *  \brief Free the memory for the density array. */
// //   void FreeMemory_CPU(void);
// //   // void FreeMemory_GPU(void);
// //
// // };
// //
//
// #endif //PARTICLES
