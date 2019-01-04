
#ifdef CPU_TIME

#ifndef TIMING_FUNCTIONS_H
#define TIMING_FUNCTIONS_H

#include "global.h"


class Time
{
public:

  int n_steps;

  // Real start_step;
  // Real end_step;
  // Real time_step;

  Real time_start;
  Real time_end;
  Real time;

  #ifdef GRAVITY_CPU
  Real time_dt_min;
  Real time_dt_max;
  Real time_dt_mean;
  Real time_dt_all;

  Real time_bound_pot_min;
  Real time_bound_pot_max;
  Real time_bound_pot_mean;
  Real time_bound_pot_all;
  #endif

  Real time_hydro_min;
  Real time_hydro_max;
  Real time_hydro_mean;
  Real time_hydro_all;

  Real time_bound_min;
  Real time_bound_max;
  Real time_bound_mean;
  Real time_bound_all;


  #ifdef GRAVITY
  Real time_potential_min;
  Real time_potential_max;
  Real time_potential_mean;
  Real time_potential_all;

  #ifdef PARTICLES

  Real time_part_dens_min;
  Real time_part_dens_max;
  Real time_part_dens_mean;
  Real time_part_dens_all;

  Real time_part_dens_transf_min;
  Real time_part_dens_transf_max;
  Real time_part_dens_transf_mean;
  Real time_part_dens_transf_all;

  Real time_part_tranf_min;
  Real time_part_tranf_max;
  Real time_part_tranf_mean;
  Real time_part_tranf_all;


  Real time_advance_particles_1_min;
  Real time_advance_particles_1_max;
  Real time_advance_particles_1_mean;
  Real time_advance_particles_1_all;

  Real time_advance_particles_2_min;
  Real time_advance_particles_2_max;
  Real time_advance_particles_2_mean;
  Real time_advance_particles_2_all;
  #endif
  #endif

  Time();
  void Initialize();
  void Start_Timer();
  void End_and_Record_Time( int time_var );
  void Print_Times();
  void Get_Average_Times();
  void Print_Average_Times();
};


#endif
#endif
