#ifdef  COSMOLOGY
#ifdef COOLING_GRACKLE

#ifndef INIT_GRACKLE_H
#define INIT_GRACKLE_H

#include"../global.h"
#include"../cosmology/cosmology.h"
#include"../gravity/grav3D.h"

extern "C" {
#include <grackle.h>
}

class Cool_GK
{
  public:


  code_units units;
  chemistry_data *data;

  Real dens_conv;
  Real momentum_conv;
  Real energy_conv;

  Real dens_to_CGS;
  Real vel_to_CGS;
  Real energy_to_CGS;

  Real gamma;

  Real temperature_units;

  #ifdef OUTPUT_TEMPERATURE
  Real *temperature;
  #endif

  #ifdef OUTPUT_COOLING_RATE
  Real *cooling_time;
  Real *cooling_rate;
  #endif

  Real tiny_number;

  // Create struct for storing grackle field data
  grackle_field_data fields;
  int field_size;




  // int procID_pfft;
  // int nproc_pfft;
  // MPI_Comm comm_pfft;


Cool_GK( void );

// void Do_Cooling_Step( Real dt );

};
void Initialize_Grackle( Cool_GK &Cool, struct parameters P,  Grav3D &Grav, Cosmology &Cosmo );

#endif
#endif
#endif
