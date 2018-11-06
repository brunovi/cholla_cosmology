#ifdef PARTICLES

#include <unistd.h>
#include "particles_3D.h"
#include "../io.h"
#include "../random_functions.h"
#include "io_particles.h"
#include "../global.h"


Particles_3D::Particles_3D( void ){}


void Particles_3D::Initialize( struct parameters P, Grav3D &Grav, Real xblocal, Real yblocal, Real zblocal, Real xbound, Real ybound, Real zbound, Real xdglobal, Real ydglobal, Real zdglobal ){

  chprintf( "\nInitializing Particles...\n");
  n_local = 0;
  n_total = 0;

  dt = 0.0;
  t = 0.0;
  dt_max = 0.01;

  C_cfl = 0.05;

  int_vector_t partIDs;
  real_vector_t mass;
  real_vector_t pos_x;
  real_vector_t pos_y;
  real_vector_t pos_z;
  real_vector_t vel_x;
  real_vector_t vel_y;
  real_vector_t vel_z;
  real_vector_t grav_x;
  real_vector_t grav_y;
  real_vector_t grav_z;

  #ifdef MPI_CHOLLA
  int_vector_t out_indxs_vec_x0;
  int_vector_t out_indxs_vec_x1;
  int_vector_t out_indxs_vec_y0;
  int_vector_t out_indxs_vec_y1;
  int_vector_t out_indxs_vec_z0;
  int_vector_t out_indxs_vec_z1;
  #endif

  G.nx_local = Grav.nx_local;
  G.ny_local = Grav.ny_local;
  G.nz_local = Grav.nz_local;
  G.nx_total = Grav.nx_total;
  G.ny_total = Grav.ny_total;
  G.nz_total = Grav.nz_total;

  G.dx = Grav.dx;
  G.dy = Grav.dy;
  G.dz = Grav.dz;

  G.xMin = xblocal;
  G.yMin = yblocal;
  G.zMin = zblocal;

  G.xMax = G.xMin + G.nx_local*G.dx;
  G.yMax = G.yMin + G.ny_local*G.dy;
  G.zMax = G.zMin + G.nz_local*G.dz;

  G.domainMin_x = xbound;
  G.domainMax_x = xbound + xdglobal;
  G.domainMin_y = ybound;
  G.domainMax_y = ybound + ydglobal;
  G.domainMin_z = zbound;
  G.domainMax_z = zbound + zdglobal;

  G.n_ghost_particles_grid = 1;
  G.n_cells = (G.nx_local+2*G.n_ghost_particles_grid) * (G.ny_local+2*G.n_ghost_particles_grid) * (G.nz_local+2*G.n_ghost_particles_grid);

  INITIAL = true;

  AllocateMemory_CPU();

  Initialize_values_CPU();

  if (strcmp(P.init, "Sphere_collapse_3D")==0) Initialize_Sphere();

  if (strcmp(P.init, "Uniform_Grid")==0) Initialize_Uniform_Grid();

  if (strcmp(P.init, "Read_Grid")==0)  Load_Particles_Data( *this, P );

  chprintf("Particles Initialized: \n n_local: %lu \n", n_local );
  chprintf(" xDomain_local:  [%.4f %.4f ] [%.4f %.4f ] [%.4f %.4f ]\n", G.xMin, G.xMax, G.yMin, G.yMax, G.zMin, G.zMax );
  chprintf(" xDomain_global: [%.4f %.4f ] [%.4f %.4f ] [%.4f %.4f ]\n", G.domainMin_x, G.domainMax_x, G.domainMin_y, G.domainMax_y, G.domainMin_z, G.domainMax_z);
  chprintf(" dx: %f  %f  %f\n", G.dx, G.dy, G.dz );

  // #ifdef MPI_CHOLLA
  // for ( int i=0; i<nproc; i++ ){
  //   if ( procID == i ) std::cout << "  [pId:"  << procID << "]  x[ "<<G.xMin<<" "<<G.xMax<<" ]" <<" y[ "<<G.yMin<<" "<<G.yMax<<" ]" <<" z[ "<<G.zMin<<" "<<G.zMax<<" ]" << std::endl;
  //   usleep(500);
  //   MPI_Barrier(world);
  // }
  // MPI_Barrier(world);
  // #endif



  #ifdef PARTICLES_OMP
  chprintf(" Using OMP for particles calculations\n");
  int n_omp_max = omp_get_max_threads();
  chprintf("  MAX OMP Threads: %d\n", n_omp_max);
  chprintf("  N OMP Threads: %d\n", N_OMP_PARTICLE_THREADS);
  // omp_set_num_threads(8);
  #pragma omp parallel num_threads( N_OMP_PARTICLE_THREADS )
  {
    int omp_id, n_omp_procs;
    part_int_t omp_pIndx_start, omp_pIndx_end;
    int omp_gridIndx_start, omp_gridIndx_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    #pragma omp barrier
    Get_OMP_Indxs( n_local, n_omp_procs, omp_id, G.nz_local + 2*G.n_ghost_particles_grid, &omp_pIndx_start, &omp_pIndx_end, &omp_gridIndx_start, &omp_gridIndx_end );

    for (int omp_indx = 0; omp_indx<n_omp_procs; omp_indx++){
      if (omp_id == omp_indx) chprintf( "  omp_id:%d  p_start:%ld  p_end:%ld  g_start:%d  g_end:%d\n", omp_id, omp_pIndx_start, omp_pIndx_end, omp_gridIndx_start, omp_gridIndx_end );
    }
  }
  #endif

  #ifdef MPI_CHOLLA
  // int n_max = std::max( G.nx_local, G.ny_local );
  // n_max = std::max( n_max, G.nz_local );
  // N_PARTICLES_TRANSFER = n_max * n_max * 10 ;
  // N_DATA_PER_PARTICLE_TRANSFER = 8;
  // N_HEADER_PARTICLES_TRANSFER = 1;
  chprintf( " N_Particles Boundaries Buffer Size: %d\n", N_PARTICLES_TRANSFER);
  // chprintf( " N_Particles Boundaries Secondary Buffer Size: %d\n", N_PARTICLES_TRANSFER_SECONDARY);
  chprintf( " N_Data per Particle Transfer: %d\n", N_DATA_PER_PARTICLE_TRANSFER);
  chprintf( " N_Header Particle Transfer: %d\n", N_HEADER_PARTICLES_TRANSFER);
  #endif

  chprintf( "\n");

}

void Particles_3D::AllocateMemory_CPU( void ){
  G.density   = (Real *) malloc(G.n_cells*sizeof(Real));
  G.gravity_x = (Real *) malloc(G.n_cells*sizeof(Real));
  G.gravity_y = (Real *) malloc(G.n_cells*sizeof(Real));
  G.gravity_z = (Real *) malloc(G.n_cells*sizeof(Real));
}

void Particles_3D::FreeMemory_CPU(void)
{
 free(G.density);
 free(G.gravity_x);
 free(G.gravity_y);
 free(G.gravity_z);

 partIDs.clear();
 pos_x.clear();
 pos_y.clear();
 pos_z.clear();
 vel_x.clear();
 vel_y.clear();
 vel_z.clear();
 grav_x.clear();
 grav_y.clear();
 grav_z.clear();

 // free(send_buffer_x0_particles);
 // free(send_buffer_x1_particles);
 // free(recv_buffer_x0_particles);
 // free(recv_buffer_x1_particles);
 // free(send_buffer_y0_particles);
 // free(send_buffer_y1_particles);
 // free(recv_buffer_y0_particles);
 // free(recv_buffer_y1_particles);
 // free(send_buffer_z0_particles);
 // free(send_buffer_z1_particles);
 // free(recv_buffer_z0_particles);
 // free(recv_buffer_z1_particles);
}

void Particles_3D::Reset( void ){
  FreeMemory_CPU();
}

void Particles_3D::Initialize_values_CPU( void ){
  int id;

  for( id=0; id<G.n_cells; id++ ){
    G.density[id] = 0;
    G.gravity_x[id] = 0;
    G.gravity_y[id] = 0;
    G.gravity_z[id] = 0;
  }
}


void Particles_3D::Initialize_Uniform_Grid( void ){
  int i, j, k;
  Real pPos_x, pPos_y, pPos_z;
  part_int_t pID = 0;
  Real Mparticle = G.dx*G.dy*G.dz;
  Real offset_x, offset_y, offset_z;
  offset_x = 0.1*G.dx;
  offset_y = -0.2*G.dy;
  offset_z = -0.3*G.dz;
  for (k=0; k<G.nz_local; k++ ){
    for (j=0; j<G.ny_local; j++ ){
      for (i=0; i<G.nx_local; i++ ){
        pPos_x = G.xMin + k*G.dx + 0.5*G.dx + offset_x;
        pPos_y = G.yMin + j*G.dy + 0.5*G.dy + offset_y;
        pPos_z = G.zMin + i*G.dz + 0.5*G.dz + offset_z;
        partIDs.push_back( pID );
        mass.push_back( Mparticle );
        pos_x.push_back( pPos_x );
        pos_y.push_back( pPos_y );
        pos_z.push_back( pPos_z);
        vel_x.push_back( 0.0 );
        vel_y.push_back( 0.0 );
        vel_z.push_back( 0.0 );
        grav_x.push_back( 0.0 );
        grav_y.push_back( 0.0 );
        grav_z.push_back( 0.0 );
        pID += 1;
      }
    }
  }
  n_local = mass.size();
  chprintf( " Particles Uniform Grid Initialized, n_local: %lu\n", n_local);
}



void Particles_3D::Initialize_Sphere( void ){
  int i, j, k, id;
  Real center_x, center_y, center_z, radius, sphereR;
  center_x = 0.5;
  center_y = 0.5;
  center_z = 0.5;;
  sphereR = 0.2;

  part_int_t n_particles_local = G.nx_local*G.ny_local*G.nz_local;
  part_int_t n_particles_total = G.nx_total*G.ny_total*G.nz_total;

  Real rho_start = 1;
  Real M_sphere = 4./3 * M_PI* rho_start * sphereR*sphereR*sphereR;
  Real Mparticle = M_sphere / n_particles_total;

  part_int_t pID = 0;
  Real pPos_x, pPos_y, pPos_z, r;
  while ( pID < n_particles_local ){
    pPos_x = Rand_Real( G.xMin, G.xMax );
    pPos_y = Rand_Real( G.yMin, G.yMax );
    pPos_z = Rand_Real( G.zMin, G.zMax );

    r = sqrt( (pPos_x-center_x)*(pPos_x-center_x) + (pPos_y-center_y)*(pPos_y-center_y) + (pPos_z-center_z)*(pPos_z-center_z) );
    if ( r > sphereR ) continue;
    partIDs.push_back( pID );
    mass.push_back( Mparticle );
    pos_x.push_back( pPos_x );
    pos_y.push_back( pPos_y );
    pos_z.push_back( pPos_z);
    vel_x.push_back( 0.0 );
    vel_y.push_back( 0.0 );
    vel_z.push_back( 0.0 );
    grav_x.push_back( 0.0 );
    grav_y.push_back( 0.0 );
    grav_z.push_back( 0.0 );
    pID += 1;
  }

  n_local = mass.size();

  chprintf( " Particles Uniform Sphere Initialized, n_local: %lu\n", n_local);

}


#endif
