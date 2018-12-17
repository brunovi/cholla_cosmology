
#ifdef CPU_TIME

#include "timing_functions.h"
#include "io.h"


#ifdef MPI_CHOLLA
#include "mpi_routines.h"
#endif

Time::Time( void ){}


void Time::Initialize(){

  time_hydro_all = 0;
  time_bound_all = 0;

  #ifdef GRAVITY
  time_potential_all = 0;

  #ifdef PARTICLES
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
    time_hydro_all += t_avg;
  }
  if( time_var == 2 ){
    time_bound_min = t_min;
    time_bound_max = t_max;
    time_bound_mean = t_avg;
    time_bound_all += t_avg;


  }
  if( time_var == 3 ){
    time_potential_min = t_min;
    time_potential_max = t_max;
    time_potential_mean = t_avg;
    time_potential_all += t_avg;
  }
}

void Time::Print_Times(){
  chprintf(" Time Hydro  min: %9.4f  max: %9.4f  avg: %9.4f\n", time_potential_min, time_potential_max, time_potential_mean);
  chprintf(" Time Boundaries  min: %9.4f  max: %9.4f  avg: %9.4f\n", time_bound_min, time_bound_max, time_bound_mean);
  #ifdef GRAVITY
  chprintf(" Time Potential  min: %9.4f  max: %9.4f  avg: %9.4f\n", time_potential_min, time_potential_max, time_potential_mean);
  #endif

}

#endif
