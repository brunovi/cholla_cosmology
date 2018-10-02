#ifdef PARTICLES

#include"density_CIC.h"
#include"../io.h"

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

void Clear_Density( Particles_3D &Parts ){
  for( int i=0; i<Parts.G.n_cells; i++ ) Parts.G.density[i] = 0;
}

void Get_Indexes_CIC( Real xMin, Real yMin, Real zMin, Real dx, Real dy, Real dz, Real pos_x, Real pos_y, Real pos_z, int &indx_x, int &indx_y, int &indx_z ){
  indx_x = (int) floor( ( pos_x - xMin - 0.5*dx ) / dx );
  indx_y = (int) floor( ( pos_y - yMin - 0.5*dy ) / dy );
  indx_z = (int) floor( ( pos_z - zMin - 0.5*dz ) / dz );
}



void Get_Density_CIC( Particles_3D &Parts ){
  int nGHST = Parts.G.n_ghost_particles_grid;
  int nx_g = Parts.G.nx_local + 2*nGHST;
  int ny_g = Parts.G.ny_local + 2*nGHST;
  int nz_g = Parts.G.nz_local + 2*nGHST;
  int n_cells = Parts.G.n_cells;

  Real xMin, yMin, zMin, dx, dy, dz;
  xMin = Parts.G.xMin;
  yMin = Parts.G.yMin;
  zMin = Parts.G.zMin;
  dx = Parts.G.dx;
  dy = Parts.G.dy;
  dz = Parts.G.dz;

  part_int_t pIndx;
  int indx_x, indx_y, indx_z, indx;
  Real pMass, x_pos, y_pos, z_pos;

  Real cell_center_x, cell_center_y, cell_center_z;
  Real delta_x, delta_y, delta_z;
  Real dV_inv = 1./(Parts.G.dx*Parts.G.dy*Parts.G.dz);
  bool ignore;
  // Real  max_pos_x =0;
  // Real dens_avrg = 0;
  for ( pIndx=0; pIndx < Parts.n_local; pIndx++ ){
    ignore = false;
    pMass = Parts.mass[pIndx] * dV_inv;
    x_pos = Parts.pos_x[pIndx];
    y_pos = Parts.pos_y[pIndx];
    z_pos = Parts.pos_z[pIndx];
    // if (x_pos > max_pos_x) max_pos_x = x_pos;
    // dens_avrg += pMass;
    Get_Indexes_CIC( xMin, yMin, zMin, dx, dy, dz, x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );
    if ( indx_x < -1 ) ignore = true;
    if ( indx_y < -1 ) ignore = true;
    if ( indx_z < -1 ) ignore = true;
    if ( indx_x > nx_g-3  ) ignore = true;
    if ( indx_y > ny_g-3  ) ignore = true;
    if ( indx_y > nz_g-3  ) ignore = true;
    if ( ignore ){
      std::cout << "ERROR CIC Index    pID: " << Parts.partIDs[pIndx] << std::endl;
      std::cout << "Negative xIndx: " << x_pos << "  " << indx_x << std::endl;
      std::cout << "Negative zIndx: " << z_pos << "  " << indx_z << std::endl;
      std::cout << "Negative yIndx: " << y_pos << "  " << indx_y << std::endl;
      std::cout << "Excess xIndx: " << x_pos << "  " << indx_x << std::endl;
      std::cout << "Excess yIndx: " << y_pos << "  " << indx_y << std::endl;
      std::cout << "Excess zIndx: " << z_pos << "  " << indx_z << std::endl;
      std::cout << std::endl;
      exit(-1);
      continue;
    }
    cell_center_x = xMin + indx_x*dx + 0.5*dx;
    cell_center_y = yMin + indx_y*dy + 0.5*dy;
    cell_center_z = zMin + indx_z*dz + 0.5*dz;
    delta_x = 1 - ( x_pos - cell_center_x ) / dx;
    delta_y = 1 - ( y_pos - cell_center_y ) / dy;
    delta_z = 1 - ( z_pos - cell_center_z ) / dz;
    indx_x += nGHST;
    indx_y += nGHST;
    indx_z += nGHST;

    indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
    Parts.G.density[indx] += pMass  * delta_x * delta_y * delta_z;

    indx = (indx_x+1) + indx_y*nx_g + indx_z*nx_g*ny_g;
    Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * delta_z;

    indx = indx_x + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
    Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * delta_z;

    indx = indx_x + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
    Parts.G.density[indx] += pMass  * delta_x * delta_y * (1-delta_z);

    indx = (indx_x+1) + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
    Parts.G.density[indx] += pMass  * (1-delta_x) * (1-delta_y) * delta_z;

    indx = (indx_x+1) + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
    Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * (1-delta_z);

    indx = indx_x + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * (1-delta_z);

    indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    Parts.G.density[indx] += pMass * (1-delta_x) * (1-delta_y) * (1-delta_z);
  }
  // dens_avrg /= (Parts.G.nx_local * Parts.G.ny_local * Parts.G.nz_local );
  // std::cout << "           Dens Avg: " << dens_avrg << std::endl;
}

#ifdef PARTICLES_OMP
void Get_Density_CIC_OMP( Particles_3D &Parts ){

  #pragma omp parallel num_threads( N_OMP_PARTICLE_THREADS )
  {
    int omp_id, n_omp_procs;
    int g_start, g_end;
    int nGHST, nx_g, ny_g, nz_g;
    nGHST = Parts.G.n_ghost_particles_grid;
    nx_g = Parts.G.nx_local + 2*nGHST;
    ny_g = Parts.G.ny_local + 2*nGHST;
    nz_g = Parts.G.nz_local + 2*nGHST;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, Parts.G.nz_local,  &g_start, &g_end );

    chprintf( "omp_id: %d  g_start: %d     g_end: %d\n", omp_id, g_start, g_end );

    // g_start -= nGHST;
    // g_end   -= nGHST;

    Real xMin, yMin, zMin, dx, dy, dz, zMin_local_p, zMax_local_p, zMin_local_g, zMax_local_g;
    dx = Parts.G.dx;
    dy = Parts.G.dy;
    dz = Parts.G.dz;
    xMin = Parts.G.xMin;
    yMin = Parts.G.yMin;
    zMin = Parts.G.zMin;

    zMin_local_g = Parts.G.zMin + g_start*dz;
    zMax_local_g = Parts.G.zMin + g_end * dz;

    zMin_local_p = Parts.G.zMin + g_start*dz - 0.5*dz ;
    zMax_local_p = Parts.G.zMin + g_end * dz + 0.5*dz ;

    int indx_x, indx_y, indx_z, indx;
    Real pMass, x_pos, y_pos, z_pos;

    Real cell_center_x, cell_center_y, cell_center_z;
    Real delta_x, delta_y, delta_z;
    Real dV_inv = 1./(Parts.G.dx*Parts.G.dy*Parts.G.dz);


    bool add_1, add_2;
    part_int_t pIndx;
    for ( pIndx=0; pIndx<Parts.n_local; pIndx++ ){
      z_pos = Parts.pos_z[pIndx];
      if ( ( z_pos < zMin_local_g) || ( z_pos >= zMax_local_g ) ) continue;
      pMass = Parts.mass[pIndx] * dV_inv;
      x_pos = Parts.pos_x[pIndx];
      y_pos = Parts.pos_y[pIndx];
      Get_Indexes_CIC( xMin, yMin, zMin, dx, dy, dz, x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );


      cell_center_x = xMin + indx_x*dx + 0.5*dx;
      cell_center_y = yMin + indx_y*dy + 0.5*dy;
      cell_center_z = zMin + indx_z*dz + 0.5*dz;
      delta_x = 1 - ( x_pos - cell_center_x ) / dx;
      delta_y = 1 - ( y_pos - cell_center_y ) / dy;
      delta_z = 1 - ( z_pos - cell_center_z ) / dz;
      indx_x += nGHST;
      indx_y += nGHST;
      indx_z += nGHST;

      add_1 = true;
      add_2 = true;

      // if (z_pos < zMin_local_g  ) add_1 = false;
      // // // if (indx_z > g_end-1 ) add_1 = false;
      // if ( z_pos > zMax_local_g ) add_2 = false;
      //
      // if ( add_1 ){
      if( z_pos < zMin_local_g + 0.5*dz ){
        indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass  * delta_x * delta_y * delta_z;

        indx = (indx_x+1) + indx_y*nx_g + indx_z*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * delta_z;

        indx = indx_x + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * delta_z;

        indx = (indx_x+1) + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass  * (1-delta_x) * (1-delta_y) * delta_z;
      }
      else{
        indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
        Parts.G.density[indx] += pMass  * delta_x * delta_y * delta_z;

        indx = (indx_x+1) + indx_y*nx_g + indx_z*nx_g*ny_g;
        Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * delta_z;

        indx = indx_x + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
        Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * delta_z;

        indx = (indx_x+1) + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
        Parts.G.density[indx] += pMass  * (1-delta_x) * (1-delta_y) * delta_z;
      }
      // if ( add_2 ){
      if ( z_pos > zMax_local_g - 0.5*dz ){
        indx = indx_x + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass  * delta_x * delta_y * (1-delta_z);

        indx = (indx_x+1) + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * (1-delta_z);

        indx = indx_x + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * (1-delta_z);

        indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
        #pragma omp atomic update
        Parts.G.density[indx] += pMass * (1-delta_x) * (1-delta_y) * (1-delta_z);
      }
      else{
        indx = indx_x + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
        Parts.G.density[indx] += pMass  * delta_x * delta_y * (1-delta_z);

        indx = (indx_x+1) + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
        Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * (1-delta_z);

        indx = indx_x + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
        Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * (1-delta_z);

        indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
        Parts.G.density[indx] += pMass * (1-delta_x) * (1-delta_y) * (1-delta_z);
      }
    }

  }
}
//
// void Get_Density_CIC_OMP( Particles_3D &Parts ){
//
//   #pragma omp parallel num_threads( N_OMP_PARTICLE_THREADS )
//   {
//     int omp_id, n_omp_procs;
//     int g_start, g_end;
//     int nGHST, nx_g, ny_g, nz_g;
//     nGHST = Parts.G.n_ghost_particles_grid;
//     nx_g = Parts.G.nx_local + 2*nGHST;
//     ny_g = Parts.G.ny_local + 2*nGHST;
//     nz_g = Parts.G.nz_local + 2*nGHST;
//
//     omp_id = omp_get_thread_num();
//     n_omp_procs = omp_get_num_threads();
//     Get_OMP_Grid_Indxs( n_omp_procs, omp_id, Parts.G.nz_local+2*nGHST,  &g_start, &g_end );
//
//     g_start -= nGHST;
//     g_end   -= nGHST;
//     Real xMin, yMin, zMin, dx, dy, dz, zMin_local, zMax_local;
//     dx = Parts.G.dx;
//     dy = Parts.G.dy;
//     dz = Parts.G.dz;
//     xMin = Parts.G.xMin;
//     yMin = Parts.G.yMin;
//     zMin = Parts.G.zMin;
//
//     zMin_local = Parts.G.zMin + g_start*dz;
//     zMax_local = Parts.G.zMin + g_end * dz;
//
//     int indx_x, indx_y, indx_z, indx;
//     Real pMass, x_pos, y_pos, z_pos;
//
//     Real cell_center_x, cell_center_y, cell_center_z;
//     Real delta_x, delta_y, delta_z;
//     Real dV_inv = 1./(Parts.G.dx*Parts.G.dy*Parts.G.dz);
//
//     int_vector_t onHold_indxs;
//
//     part_int_t pIndx;
//     for ( pIndx=0; pIndx<Parts.n_local; pIndx++ ){
//       z_pos = Parts.pos_z[pIndx];
//       if ( ( z_pos < zMin_local ) || ( z_pos >= zMax_local ) ) continue;
//       pMass = Parts.mass[pIndx] * dV_inv;
//       x_pos = Parts.pos_x[pIndx];
//       y_pos = Parts.pos_y[pIndx];
//       Get_Indexes_CIC( xMin, yMin, zMin, dx, dy, dz, x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );
//       if ( (indx_z == g_end-1 ) || (indx_z == g_end-2 ) ){
//         onHold_indxs.push_back( pIndx);
//         continue;
//       }
//
//       cell_center_x = xMin + indx_x*dx + 0.5*dx;
//       cell_center_y = yMin + indx_y*dy + 0.5*dy;
//       cell_center_z = zMin + indx_z*dz + 0.5*dz;
//       delta_x = 1 - ( x_pos - cell_center_x ) / dx;
//       delta_y = 1 - ( y_pos - cell_center_y ) / dy;
//       delta_z = 1 - ( z_pos - cell_center_z ) / dz;
//       indx_x += nGHST;
//       indx_y += nGHST;
//       indx_z += nGHST;
//
//       indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * delta_y * delta_z;
//
//       indx = (indx_x+1) + indx_y*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * delta_z;
//
//       indx = indx_x + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * delta_z;
//
//       indx = indx_x + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * delta_y * (1-delta_z);
//
//       indx = (indx_x+1) + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * (1-delta_x) * (1-delta_y) * delta_z;
//
//       indx = (indx_x+1) + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * (1-delta_z);
//
//       indx = indx_x + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * (1-delta_z);
//
//       indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass * (1-delta_x) * (1-delta_y) * (1-delta_z);
//     }
//
//     part_int_t indx_p;
//     for (indx_p=0; indx_p<onHold_indxs.size(); indx_p++ ){
//       pIndx = onHold_indxs[indx_p];
//       pMass = Parts.mass[pIndx] * dV_inv;
//       x_pos = Parts.pos_x[pIndx];
//       y_pos = Parts.pos_y[pIndx];
//       z_pos = Parts.pos_z[pIndx];
//       Get_Indexes_CIC( xMin, yMin, zMin, dx, dy, dz, x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );
//
//       cell_center_x = xMin + indx_x*dx + 0.5*dx;
//       cell_center_y = yMin + indx_y*dy + 0.5*dy;
//       cell_center_z = zMin + indx_z*dz + 0.5*dz;
//       delta_x = 1 - ( x_pos - cell_center_x ) / dx;
//       delta_y = 1 - ( y_pos - cell_center_y ) / dy;
//       delta_z = 1 - ( z_pos - cell_center_z ) / dz;
//       indx_x += nGHST;
//       indx_y += nGHST;
//       indx_z += nGHST;
//
//       indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * delta_y * delta_z;
//
//       indx = (indx_x+1) + indx_y*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * delta_z;
//
//       indx = indx_x + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * delta_z;
//
//       indx = indx_x + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * delta_y * (1-delta_z);
//
//       indx = (indx_x+1) + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * (1-delta_x) * (1-delta_y) * delta_z;
//
//       indx = (indx_x+1) + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * (1-delta_x) * delta_y * (1-delta_z);
//
//       indx = indx_x + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass  * delta_x * (1-delta_y) * (1-delta_z);
//
//       indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
//       Parts.G.density[indx] += pMass * (1-delta_x) * (1-delta_y) * (1-delta_z);
//
//     }
//
//
//   }
// }

#endif

void Copy_Particles_Density_to_Gravity( Grid3D &G ){
  int nx, ny, nz, nGHST;
  nx = G.Particles.G.nx_local;
  ny = G.Particles.G.ny_local;
  nz = G.Particles.G.nz_local;
  nGHST = G.Particles.G.n_ghost_particles_grid;
  int i, j, k, id_CIC, id_grid;
  for ( k=0; k<nz; k++ ){
    for ( j=0; j<ny; j++ ){
      for ( i=0; i<nx; i++ ){
      id_CIC = (i+nGHST) + (j+nGHST)*(nx+2*nGHST) + (k+nGHST)*(nx+2*nGHST)*(ny+2*nGHST);
      id_grid = i + j*nx + k*nx*ny;
      // #ifndef COSMOLOGY
      // // G.Grav.F.density_h[id_grid] += 4 * M_PI * G.Particles.G.density[id_CIC];
      // G.Grav.F.density_h[id_grid] += G.Particles.G.density[id_CIC];
      // #else
      G.Grav.F.density_h[id_grid] += G.Particles.G.density[id_CIC];
      // G.Grav.F.density_h[id_grid] += 4 * M_PI * G.Cosmo.cosmo_G * ( G.Particles.G.density[id_CIC] - G.Cosmo.rho_0_dm ) / G.Cosmo.current_a ;
      // G.Grav.F.density_h[id_grid] += 4 * M_PI * G.Cosmo.cosmo_G * ( G.Particles.G.density[id_CIC] - G.Cosmo.dens_avrg ) / G.Cosmo.current_a ;
      // #endif
      }
    }
  }
}



void Get_Particles_Density_CIC( Grid3D &G, struct parameters P, Real *time_pDens, Real *time_pDens_trans ){

  // chprintf(" Computing Particles Density: CIC\n");

  Real start, stop, time_dens, time_dens_trans;

  Clear_Density( G.Particles );

  start = get_time();
  #ifdef PARTICLES_OMP
  // Get_Density_CIC_OMP( G.Particles );
  Get_Density_CIC( G.Particles );
  #elif PARTICLES_CUDA
  Get_Density_CIC_CUDA( G.Particles );
  #else
  Get_Density_CIC( G.Particles );
  #endif
  stop = get_time();
  time_dens = (stop - start) * 1000.0;

  start = get_time();
  #ifndef MPI_CHOLLA
  Transfer_Particles_Density_Boundaries( G.Particles );
  #else
  Transfer_Particles_Density_Boundaries_MPI( G, P );
  #endif
  stop = get_time();
  time_dens_trans = (stop - start) * 1000.0;



  *time_pDens = time_dens;
  *time_pDens_trans = time_dens_trans;

}

#endif
