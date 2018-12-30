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
#ifdef POTENTIAL_PFFT
#include "gravity/potential_PFFT_3D.h"
#endif

#ifdef PARTICLES
#include "particles/particles_dynamics.h"
#endif

#ifdef COSMOLOGY
#include "cosmology/cosmology.h"
#include "cosmology/cosmology_units.h"
#include "cosmology/io_cosmology.h"
#endif

#endif

#ifdef COOLING_GRACKLE
#include "cooling/cool_grackle.h"
#include "cooling/grackle_functions.h"
#endif


#ifdef CPU_TIME
#include "timing_functions.h"
#endif

#define OUTPUT
// #define CPU_TIME

int main(int argc, char *argv[])
{
  // timing variables
  double start_total, stop_total, start_step, stop_step;
  #ifdef CPU_TIME
  double stop_init, init_min, init_max, init_avg;
  // double start_bound, stop_bound, bound_min, bound_max, bound_avg;
  // double start_hydro, stop_hydro, hydro_min, hydro_max, hydro_avg;
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
  G.Grav.Initialize( G.H.xblocal, G.H.yblocal, G.H.zblocal, G.H.xdglobal, G.H.ydglobal, G.H.zdglobal, P.nx, P.ny, P.nz, G.H.nx_real, G.H.ny_real, G.H.nz_real, G.H.dx, G.H.dy, G.H.dz, G.H.n_ghost_pot_offset  );

  #ifdef POTENTIAL_CUFFT
  Potential_CUFFT_3D p_solver;
  #endif
  #ifdef POTENTIAL_FFTW
  Potential_FFTW_3D p_solver;
  #endif
  #ifdef POTENTIAL_PFFT
  Potential_PFFT_3D p_solver;
  #endif
  p_solver.Initialize( G.Grav );
  #endif

  #ifdef PARTICLES
  G.Particles.Initialize( P, G.Grav, G.H.xblocal, G.H.yblocal, G.H.zblocal, G.H.xbound, G.H.ybound, G.H.zbound, G.H.xdglobal, G.H.ydglobal, G.H.zdglobal );
  #endif

  #ifdef COSMOLOGY
  Initialize_Cosmology( G.Cosmo, P, G.Particles, G.Grav );
  Change_Cosmological_Frame_Sytem( G, true );
  #endif

  #ifdef GRAVITY
  Compute_Gravitational_Potential( G, p_solver,  P );
  Copy_Potential_To_Hydro_Grid( G );
  #endif


  // set boundary conditions (assign appropriate values to ghost cells)
  chprintf("\nSetting boundary conditions...\n");
  G.H.TRANSFER_HYDRO_BOUNDARIES = true;
  G.Set_Boundary_Conditions(P);
  G.H.TRANSFER_HYDRO_BOUNDARIES = false;
  chprintf("Boundary conditions set.\n");

  #ifdef PARTICLES
  //Transfer Particles Boundaries
  G.Particles.TRANSFER_PARTICLES_BOUNDARIES = true;
  G.Set_Boundary_Conditions(P);
  G.Particles.TRANSFER_PARTICLES_BOUNDARIES = false;
  #endif

  #ifdef COOLING_GRACKLE
  Initialize_Grackle( G.Cool, P, G.Grav, G.Cosmo );
  Initialize_Grackle_Fields( G );
  Real start_cool, stop_cool, time_cool;
  start_cool = get_time();
  Do_Cooling_Step( 0, G );
  stop_cool = get_time();
  time_cool = stop_cool - start_cool;
  chprintf( " Time Cooling: %f\n", time_cool*1000 );
  #endif

  #ifdef PARTICLES
  Get_Particles_Acceleration( G, 0, G.Particles.n_local, 0, G.Particles.G.nz_local + 2*G.Particles.G.n_ghost_particles_grid );
  #endif


  chprintf("Dimensions of each cell: dx = %f dy = %f dz = %f\n", G.H.dx, G.H.dy, G.H.dz);
  chprintf("Ratio of specific heats gamma = %f\n",gama);
  chprintf("Nstep = %d  Timestep = %f  Simulation time = %f\n", G.H.n_step, G.H.dt, G.H.t);


  #ifdef OUTPUT
  // write the initial conditions to file
  chprintf("Writing initial conditions to file...\n");
  WriteData(G, P, nfile);
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

  #ifdef CPU_TIME
  G.Timer.Initialize();
  #endif

  bool output_now = false;
  // Evolve the grid, one timestep at a time
  chprintf("Starting calculations.\n");
  while (G.H.t < P.tout)
  //while (G.H.n_step < 1)
  {
    #ifdef COOLING_GRACKLE
    break;
    #endif

    chprintf("n_step: %d \n", G.H.n_step + 1 );
    // get the start time
    start_step = get_time();

    // calculate the timestep
    G.set_dt(dti);

    #ifndef COSMOLOGY
    if (G.H.t + G.H.dt > outtime)
    {
      G.H.dt = outtime - G.H.t;
    }

    #ifdef MPI_CHOLLA
    G.H.dt = ReduceRealMin(G.H.dt);
    #endif
    #endif

    #ifdef GRAVITY
    Set_dt( G, output_now, G.H.n_step + 1 );
    // Extrapolate gravitational potential for hydro step
    Extrapolate_Grav_Potential( G );
    #endif


    // Advance the grid by one timestep
    #ifdef CPU_TIME
    G.Timer.Start_Timer();
    dti = G.Update_Grid();
    G.Timer.End_and_Record_Time( 1 );
    #else
    dti = G.Update_Grid();
    #endif //CPU_TIME


    #ifdef PARTICLES
    //Advance the particles ( first step )
    Update_Particles( G, 1 );
    #endif

    // update the simulation time ( t += dt )
    G.Update_Time();

    // add one to the timestep count
    G.H.n_step++;

    //Compute Gravitational potential for next step
    #ifdef GRAVITY
    Compute_Gravitational_Potential( G, p_solver, P );
    Copy_Potential_To_Hydro_Grid( G );
    #endif

    #ifdef GRAVITY_CORRECTOR
    Apply_Gavity_Corrector( G, P);
    #endif

    #ifdef REVERT_STEP
    Get_Delta_Conserved( G );
    #endif

    // set boundary conditions for next time step
    #ifdef CPU_TIME
    G.Timer.Start_Timer();
    G.H.TRANSFER_HYDRO_BOUNDARIES = true;
    G.Set_Boundary_Conditions(P);
    G.H.TRANSFER_HYDRO_BOUNDARIES = false;
    G.Timer.End_and_Record_Time( 2 );
    #else
    G.Set_Boundary_Conditions(P);
    #endif //CPU_TIME

    #ifdef PARTICLES
    //Transfer Particles Boundaries
    G.Particles.TRANSFER_PARTICLES_BOUNDARIES = true;
    G.Timer.Start_Timer();
    G.Set_Boundary_Conditions(P);
    G.Timer.End_and_Record_Time( 8 );
    G.Particles.TRANSFER_PARTICLES_BOUNDARIES = false;
    #endif

    #ifdef PARTICLES
    //Advance the particles ( second step )
    Update_Particles( G, 2 );
    #endif

    // #ifdef CPU_TIME
    // #ifdef MPI_CHOLLA
    // chprintf("hydro min: %9.4f  max: %9.4f  avg: %9.4f\n", hydro_min, hydro_max, hydro_avg);
    // chprintf("bound min: %9.4f  max: %9.4f  avg: %9.4f\n", bound_min, bound_max, bound_avg);
    // #endif //MPI_CHOLLA
    // #endif //CPU_TIME

    #ifdef CPU_TIME
    G.Timer.Print_Times();
    #endif

    // get the time to compute the total timestep
    stop_step = get_time();
    stop_total = get_time();
    G.H.t_wall = stop_total-start_total;
    #ifdef MPI_CHOLLA
    G.H.t_wall = ReduceRealMax(G.H.t_wall);
    #endif
    chprintf("n_step: %d   sim time: %10.7f   sim timestep: %7.4e  timestep time = %9.3f ms   total time = %9.4f s\n\n",
      G.H.n_step, G.H.t, G.H.dt, (stop_step-start_step)*1000, G.H.t_wall);

    #ifndef COSMOLOGY
    if (G.H.t == outtime) output_now = true;
    #endif

    if ( output_now ){
      #ifdef OUTPUT
      /*output the grid data*/
      WriteData(G, P, nfile);
      // add one to the output file count
      nfile++;
      #endif //OUTPUT
      // update to the next output time
      outtime += P.outstep;

      output_now = false;

    }


    // if (G.H.n_step == 1) break;

    #ifdef COSMOLOGY
    if ( G.Cosmo.current_a >= G.Cosmo.scale_outputs[G.Cosmo.n_outputs-1] ) {
      chprintf( "\nReached Last Cosmological Output: Ending Simulation\n");
      break;
    }
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

  //
  // #ifdef CPU_TIME
  // time_hydro_total /= ( step_counter -1 );
  // time_boundaries_total /= ( step_counter -1 );
  // time_all_total /= ( step_counter -1 );
  // time_potential_total /= ( step_counter -1 );
  // time_particles_total /= ( step_counter -1 );
  // chprintf("\n\nSimulation Finished\n");
  // chprintf(" N Steps: %d\n", step_counter);
  // chprintf(" Time Average Hydro: %f\n", time_hydro_total);
  // chprintf(" Time Average Boundaries: %f\n", time_boundaries_total);
  // chprintf(" Time Average Potential: %f\n", time_potential_total);
  // chprintf(" Time Average Particles: %f\n", time_particles_total);
  // chprintf(" Time Average Total: %f\n", time_all_total);
  //
  //
  // // Output timing values
  // ofstream out_file;
  // out_file.open("run_timing.log", ios::app);
  // out_file << P.nz << " " << P.ny << " " << P.nx << " ";
  // #ifndef GRAVITY_OMP
  // out_file << 0 << " ";
  // #endif
  // #ifdef GRAVITY_OMP
  // out_file << N_OMP_GRAVITY_THREADS << " ";
  // #endif
  // out_file << time_hydro_total << " " << time_boundaries_total << " ";
  // out_file << time_potential_total << " " << time_particles_total << " ";
  // out_file << "\n";
  // out_file.close();
  // #endif


  // free the grid
  G.Reset();

  #ifdef GRAVITY
  p_solver.Reset();
  #endif

  #ifdef Particles_3D
  G.Particles.Reset();
  #endif

  #ifdef COOLING_GRACKLE
  Clear_Data_Grackle( G );
  #endif

  #ifdef MPI_CHOLLA
  MPI_Finalize();
  #endif /*MPI_CHOLLA*/

  return 0;

}
