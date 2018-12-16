
#ifdef CPU_TIME

#include "timing_functions.h"


void Time::Start_Timer(){
  time_start = get_time();
}

void Time::End_and_Record_Time( int time_var ){

  time_end = get_time();
  time = (time_end - time_start)*1000;

  if( time_var == 1 ){

  }


}
#endif
