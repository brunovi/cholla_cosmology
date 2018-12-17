
#ifdef CPU_TIME

#include "timing_functions.h"
#include "io.h"


#ifdef MPI_CHOLLA
#include "mpi_routines.h"
#endif

Time::Time( void ){}


void Time::Initialize(){

  n_steps = 0;

  time_hydro_all = 0;
  time_bound_all = 0;

  #ifdef GRAVITY
  time_potential_all = 0;

  #ifdef PARTICLES
  time_part_dens_transf_all = 0;
  time_part_dens_transf_all = 0;
  time_advance_particles_1_all = 0;
  time_advance_particles_2_all = 0;
  #endif
  #endif

  chprintf( "\nTiming Functions is ON \n\n");


}

void Time::Start_Timer(){
  time_start = get_time();
}

void Time::End_and_Record_Time( int time_var ){

  time_end = get_time();
  time = (time_end - time_start)*1000;

  Real t_min, t_max, t_avg;

  #ifdef MPI_CHOLLA
  t_min = ReduceRealMin(time);
  t_max = ReduceRealMax(time);
  t_avg = ReduceRealAvg(time);
  #endif


  if( time_var == 1 ){
    time_hydro_min = t_min;
    time_hydro_max = t_max;
    time_hydro_mean = t_avg;
    if (n_steps > 0) time_hydro_all += t_avg;
  }
  if( time_var == 2 ){
    time_bound_min = t_min;
    time_bound_max = t_max;
    time_bound_mean = t_avg;
    if (n_steps > 0) time_bound_all += t_avg;
  }
  if( time_var == 3 ){
    time_potential_min = t_min;
    time_potential_max = t_max;
    time_potential_mean = t_avg;
    if (n_steps > 0) time_potential_all += t_avg;
  }

  if( time_var == 4 ){
    time_part_dens_min = t_min;
    time_part_dens_max = t_max;
    time_part_dens_mean = t_avg;
    if (n_steps > 0) time_part_dens_all += t_avg;
  }

  if( time_var == 5 ){
    time_part_dens_transf_min = t_min;
    time_part_dens_transf_max = t_max;
    time_part_dens_transf_mean = t_avg;
    if (n_steps > 0) time_part_dens_transf_all += t_avg;
  }

  if( time_var == 6 ){
    time_advance_particles_1_min = t_min;
    time_advance_particles_1_max = t_max;
    time_advance_particles_1_mean = t_avg;
    if (n_steps > 0) time_advance_particles_1_all += t_avg;
  }

  if( time_var == 7 ){
    time_advance_particles_2_min = t_min;
    time_advance_particles_2_max = t_max;
    time_advance_particles_2_mean = t_avg;
    if (n_steps > 0) time_advance_particles_2_all += t_avg;
  }


  if ( time_var == 1 ) n_steps += 1;
}

void Time::Print_Times(){
  chprintf(" Time Hydro             min: %9.4f  max: %9.4f  avg: %9.4f\n", time_hydro_min, time_hydro_max, time_hydro_mean);
  chprintf(" Time Boundaries        min: %9.4f  max: %9.4f  avg: %9.4f\n", time_bound_min, time_bound_max, time_bound_mean);
  #ifdef GRAVITY
  chprintf(" Time Potential         min: %9.4f  max: %9.4f  avg: %9.4f\n", time_potential_min, time_potential_max, time_potential_mean);
  #ifdef PARTICLES
  chprintf(" Time Part Density      min: %9.4f  max: %9.4f  avg: %9.4f\n", time_part_dens_min, time_part_dens_max, time_part_dens_mean);
  chprintf(" Time Part Dens Transf  min: %9.4f  max: %9.4f  avg: %9.4f\n", time_part_dens_transf_min, time_part_dens_transf_max, time_part_dens_transf_mean);
  chprintf(" Time Advance Part 1    min: %9.4f  max: %9.4f  avg: %9.4f\n", time_advance_particles_1_min, time_advance_particles_1_max, time_advance_particles_1_mean);
  chprintf(" Time Advance Part 2    min: %9.4f  max: %9.4f  avg: %9.4f\n", time_advance_particles_2_min, time_advance_particles_2_max, time_advance_particles_2_mean);
  #endif
  #endif

}

#endif
