
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

  #endif

  chprintf( "\nTiming Functions is ON \n\n");


}

void Time::Start_Timer(){
  time_start = get_time();
}

void Time::End_and_Record_Time( int time_var ){

  time_end = get_time();
  time = (time_end - time_start)*1000;

  Real t_min, t_max, t_mean;

  #ifdef MPI_CHOLLA

  #endif

  if( time_var == 1 ){

  }


}
#endif
