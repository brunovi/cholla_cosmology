
#ifdef CPU_TIME

#include "global.h"

class Time
{
public:

  Real time_start;
  Real time_end;
  Real time;

  Real time_hydro_min;
  Real time_hydro_max;
  Real time_hydro_mean;

  Real time_bound_min;
  Real time_bound_max;
  Real time_bound_mean;


  #ifdef GRAVITY
  Real time_potential_min;
  Real time_potential_max;
  Real time_potential_mean;
  #endif


  void Start_Timer();
  void End_and_Record_Time( int time_var );
};


#endif
