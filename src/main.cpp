/*! \file main.cpp
 *  \brief Program to run the grid code. */

#ifdef MPI_CHOLLA
#include <mpi.h>
#include "mpi_routines.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "global.h"
#include "grid3D.h"
#include "io.h"
#include "error_handling.h"

#include <iostream>
#include <fstream>
using namespace std;

#ifdef GRAVITY
#include "gravity/gravity_functions.h"
#ifdef POTENTIAL_CUFFT
#include "gravity/potential_CUFFT_3D.h"
#endif
#ifdef POTENTIAL_FFTW
#include "gravity/potential_FFTW_3D.h"
#endif

#ifdef PARTICLES
#include "particles/particles_dynamics.h"
#endif

#endif

#define OUTPUT
#define CPU_TIME

int main(int argc, char *argv[])
{
  // timing variables
  double start_total, stop_total, start_step, stop_step;
  #ifdef CPU_TIME
  double stop_init, init_min, init_max, init_avg;
  double start_bound, stop_bound, bound_min, bound_max, bound_avg;
  double start_hydro, stop_hydro, hydro_min, hydro_max, hydro_avg;
  double init, bound, hydro;
  init = bound = hydro = 0;
  #endif //CPU_TIME

  // start the total time
  start_total = get_time();

  /* Initialize MPI communication */
  #ifdef MPI_CHOLLA
  InitializeChollaMPI(&argc, &argv);
  #endif /*MPI_CHOLLA*/

  Real dti = 0; // inverse time step, 1.0 / dt

  // input parameter variables
  char *param_file;
  struct parameters P;
  int nfile = 0; // number of output files
  Real outtime = 0; // current output time


  // read in command line arguments
  if (argc != 2)
  {
    chprintf("usage: %s <parameter_file>\n", argv[0]);
    chexit(-1);
  } else {
    param_file = argv[1];
  }

  // create the grid
  Grid3D G;

  // read in the parameters
  parse_params (param_file, &P);
  // and output to screen
  chprintf ("Parameter values:  nx = %d, ny = %d, nz = %d, tout = %f, init = %s, boundaries = %d %d %d %d %d %d\n",
    P.nx, P.ny, P.nz, P.tout, P.init, P.xl_bcnd, P.xu_bcnd, P.yl_bcnd, P.yu_bcnd, P.zl_bcnd, P.zu_bcnd);
  chprintf ("Output directory:  %s\n", P.outdir);


  // initialize the grid
  G.Initialize(&P);
  chprintf("Local number of grid cells: %d %d %d %d\n", G.H.nx_real, G.H.ny_real, G.H.nz_real, G.H.n_cells);


  // Set initial conditions and calculate first dt
  chprintf("Setting initial conditions...\n");
  G.Set_Initial_Conditions(P);
  chprintf("Initial conditions set.\n");
  // set main variables for Read_Grid inital conditions
  if (strcmp(P.init, "Read_Grid") == 0) {
    dti = C_cfl / G.H.dt;
    outtime += G.H.t;
    nfile = P.nfile*P.nfull;
  }


  #ifdef GRAVITY
  G.Grav.Initialize( G.H.xdglobal, G.H.ydglobal, G.H.zdglobal, P.nx, P.ny, P.nz, G.H.nx_real, G.H.ny_real, G.H.nz_real, G.H.dx, G.H.dy, G.H.dz, G.H.n_ghost_pot_offset  );

  #ifdef POTENTIAL_CUFFT
  Potential_CUFFT_3D p_solver;
  #endif
  #ifdef POTENTIAL_FFTW
  Potential_FFTW_3D p_solver;
  #endif

  p_solver.Initialize( G.Grav );
  #endif


  Real time_advance_particles = 0 ;
  #ifdef PARTICLES
  G.Particles.Initialize( G.Grav, G.H.xblocal, G.H.yblocal, G.H.zblocal, G.H.xbound, G.H.ybound, G.H.zbound, G.H.xdglobal, G.H.ydglobal, G.H.zdglobal );
  #endif



  Real time_potential, time_particles_density;
  time_potential = time_particles_density = 0;
  #ifdef GRAVITY
  Compute_Gravitational_Potential( G, p_solver, &time_potential, &time_particles_density );
  Copy_Potential_To_Hydro_Grid( G );
  #endif


  // set boundary conditions (assign appropriate values to ghost cells)
  chprintf("Setting boundary conditions...\n");
  G.Set_Boundary_Conditions(P);
  chprintf("Boundary conditions set.\n");

  #ifdef PARTICLES
  Get_Particles_Acceleration( G, 0, G.Particles.n_local, 0, G.Particles.G.nz_local + 2*G.Particles.G.n_ghost_particles_grid );
  #endif

  chprintf("Dimensions of each cell: dx = %f dy = %f dz = %f\n", G.H.dx, G.H.dy, G.H.dz);
  chprintf("Ratio of specific heats gamma = %f\n",gama);
  chprintf("Nstep = %d  Timestep = %f  Simulation time = %f\n", G.H.n_step, G.H.dt, G.H.t);


  #ifdef OUTPUT
  if (strcmp(P.init, "Read_Grid") != 0) {
  // write the initial conditions to file
  chprintf("Writing initial conditions to file...\n");
  WriteData(G, P, nfile);
  }
  // add one to the output file count
  nfile++;
  #endif //OUTPUT
  // increment the next output time
  outtime += P.outstep;

  #ifdef CPU_TIME
  stop_init = get_time();
  init = stop_init - start_total;
  #ifdef MPI_CHOLLA
  init_min = ReduceRealMin(init);
  init_max = ReduceRealMax(init);
  init_avg = ReduceRealAvg(init);
  chprintf("Init  min: %9.4f  max: %9.4f  avg: %9.4f\n", init_min, init_max, init_avg);
  #else
  printf("Init %9.4f\n", init);
  #endif //MPI_CHOLLA
  #endif //CPU_TIME

  #ifdef PARTICLES
  Real dt_particles, dt_min;
  #endif

  #ifdef CPU_TIME
  int step_counter = 0;
  Real time_hydro_total, time_potential_total, time_particles_total, time_all_total, time_boundaries_total;
  time_hydro_total = time_potential_total = time_particles_total = time_all_total, time_boundaries_total = 0;
  #endif

  // Evolve the grid, one timestep at a time
  chprintf("Starting calculations.\n");
  while (G.H.t < P.tout)
  //while (G.H.n_step < 1)
  {

    // get the start time
    start_step = get_time();

    // calculate the timestep
    G.set_dt(dti);

    if (G.H.t + G.H.dt > outtime)
    {
      G.H.dt = outtime - G.H.t;
    }

    #ifdef MPI_CHOLLA
    G.H.dt = ReduceRealMin(G.H.dt);
    #endif

    #ifdef PARTICLES
    dt_particles = Get_Particles_dt( G.Particles );
    dt_min = std::min( G.H.dt, dt_particles );
    G.H.dt = dt_min;
    G.Particles.dt = dt_min;
    // G.Particles.dt = 0;
    #endif

    #ifdef GRAVITY
    if ( G.Grav.INITIAL ){
      G.Grav.dt_prev = G.H.dt;
      G.Grav.dt_now = G.H.dt;
      // G.Grav.INITIAL = false;
    }else{
      G.Grav.dt_prev = G.Grav.dt_now;
      G.Grav.dt_now = G.H.dt;
    }
    #endif

    #ifdef PARTICLES
    //Advance the particles ( first step )
    time_advance_particles = Update_Particles( G, 1 );
    #endif

    // Extrapolate gravitational potential for hydro step
    #ifdef GRAVITY
    Extrapolate_Grav_Potential( G );
    #endif

    // Advance the grid by one timestep
    #ifdef CPU_TIME
    start_hydro = get_time();
    #endif //CPU_TIME
    dti = G.Update_Grid();
    #ifdef CPU_TIME
    stop_hydro = get_time();
    hydro = stop_hydro - start_hydro;
    chprintf( " Time Hydro: %f\n", hydro*1000 );
    #ifdef MPI_CHOLLA
    hydro_min = ReduceRealMin(hydro);
    hydro_max = ReduceRealMax(hydro);
    hydro_avg = ReduceRealAvg(hydro);
    #endif //MPI_CHOLLA
    #endif //CPU_TIME
    //printf("%d After Grid Update: %f %f\n", procID, G.H.dt, dti);

    // update the time
    G.H.t += G.H.dt;

    // add one to the timestep count
    G.H.n_step++;

    //Compute Gravitational potential for next step
    #ifdef GRAVITY
    Compute_Gravitational_Potential( G, p_solver, &time_potential, &time_particles_density );
    Copy_Potential_To_Hydro_Grid( G );
    #ifdef CPU_TIME
    chprintf( " Time Potential: %f\n", time_potential );
    #endif
    #endif

    // set boundary conditions for next time step
    #ifdef CPU_TIME
    start_bound = get_time();
    #endif //CPU_TIME
    G.Set_Boundary_Conditions(P);
    #ifdef CPU_TIME
    stop_bound = get_time();
    bound = stop_bound - start_bound;
    chprintf( " Time Boundaries: %f\n", bound*1000 );
    #ifdef MPI_CHOLLA
    bound_min = ReduceRealMin(bound);
    bound_max = ReduceRealMax(bound);
    bound_avg = ReduceRealAvg(bound);
    #endif //MPI_CHOLLA
    #endif //CPU_TIME

    #ifdef PARTICLES
    //Advance the particles ( second step )
    time_advance_particles += Update_Particles( G, 2 );
    #ifdef CPU_TIME
    chprintf( " Time Particles Density: %f\n", time_particles_density );
    chprintf( " Time Advance Particles: %f\n", time_advance_particles );
    chprintf( " N Local Particles: %ld\n", G.Particles.n_local );
    #endif
    #endif

    #ifdef CPU_TIME
    #ifdef MPI_CHOLLA
    chprintf("hydro min: %9.4f  max: %9.4f  avg: %9.4f\n", hydro_min, hydro_max, hydro_avg);
    chprintf("bound min: %9.4f  max: %9.4f  avg: %9.4f\n", bound_min, bound_max, bound_avg);
    #endif //MPI_CHOLLA
    #endif //CPU_TIME


    // get the time to compute the total timestep
    stop_step = get_time();
    stop_total = get_time();
    G.H.t_wall = stop_total-start_total;
    #ifdef MPI_CHOLLA
    G.H.t_wall = ReduceRealMax(G.H.t_wall);
    #endif
    chprintf("n_step: %d   sim time: %10.7f   sim timestep: %7.4e  timestep time = %9.3f ms   total time = %9.4f s\n",
      G.H.n_step, G.H.t, G.H.dt, (stop_step-start_step)*1000, G.H.t_wall);

    if (G.H.t == outtime)
    {
      #ifdef OUTPUT
      /*output the grid data*/
      WriteData(G, P, nfile);
      // add one to the output file count
      nfile++;
      #endif //OUTPUT
      // update to the next output time
      outtime += P.outstep;
    }

    #ifdef CPU_TIME
    if (step_counter > 0 ){
      time_potential_total += time_potential;
      time_particles_total += time_particles_density;
      time_particles_total += time_advance_particles;
      time_hydro_total += hydro*1000;
      time_boundaries_total += bound*1000;
      time_all_total += (stop_step-start_step)*1000;
    }
    step_counter += 1;
    #endif
/*
    // check for failures
    for (int i=G.H.n_ghost; i<G.H.nx-G.H.n_ghost; i++) {
      for (int j=G.H.n_ghost; j<G.H.ny-G.H.n_ghost; j++) {
        for (int k=G.H.n_ghost; k<G.H.nz-G.H.n_ghost; k++) {
          int id = i + j*G.H.nx + k*G.H.nx*G.H.ny;
          if (G.C.density[id] < 0.0 || G.C.density[id] != G.C.density[id]) {
            printf("Failure in cell %d %d %d. Density %e\n", i, j, k, G.C.density[id]);
            #ifdef MPI_CHOLLA
            MPI_Finalize();
            chexit(-1);
            #endif
            exit(0);
          }
        }
      }
    }
*/

  } /*end loop over timesteps*/


  #ifdef CPU_TIME
  time_hydro_total /= ( step_counter -1 );
  time_boundaries_total /= ( step_counter -1 );
  time_all_total /= ( step_counter -1 );
  time_potential_total /= ( step_counter -1 );
  time_particles_total /= ( step_counter -1 );
  chprintf("\n\nSimulation Finished\n");
  chprintf(" N Steps: %d\n", step_counter);
  chprintf(" Time Average Hydro: %f\n", time_hydro_total);
  chprintf(" Time Average Boundaries: %f\n", time_boundaries_total);
  chprintf(" Time Average Potential: %f\n", time_potential_total);
  chprintf(" Time Average Particles: %f\n", time_particles_total);
  chprintf(" Time Average Total: %f\n", time_all_total);

  // Output timing values
  ofstream out_file;
  out_file.open("run_timing.log", ios::app);
  out_file << G.Grav.nz_total << " " << G.Grav.ny_total << " " << G.Grav.nx_total << " ";
  #ifndef PARTICLES_OMP
  out_file << 0 << " ";
  #endif
  #ifdef PARTICLES_OMP
  out_file << N_OMP_PARTICLE_THREADS << " ";
  #endif
  out_file << time_hydro_total << " " << time_boundaries_total << " ";
  out_file << time_potential_total << " " << time_particles_total << " ";
  out_file << "\n";
  out_file.close();
  #endif


  // free the grid
  G.Reset();

  #ifdef GRAVITY
  p_solver.Reset();
  #endif

  #ifdef Particles_3D
  G.Particles.Reset();
  #endif

  #ifdef MPI_CHOLLA
  MPI_Finalize();
  #endif /*MPI_CHOLLA*/

  return 0;

}
