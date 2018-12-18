#ifdef PARTICLES

#ifndef PARTICLES_BOUNDARIES_H
#define PARTICLES_BOUNDARIES_H


#include "particles_3D.h"

int real_to_int( Real inVal );

void Tranfer_Particles_Boundaries( Particles_3D &Parts );

#ifdef MPI_CHOLLA

#endif



#endif
#endif
