#ifndef IO_COSMOLOGY_H
#define IO_COSMOLOGY_H

#include"../global.h"
#include"../grid3D.h"
#include"../io.h"
#include"cosmology.h"


void Load_Scale_Outputs( struct parameters P, Cosmology &Cosmo );

void Set_Next_Scale_Output( Cosmology &Cosmo );

#endif /*IO_PARTICLES_H*/
