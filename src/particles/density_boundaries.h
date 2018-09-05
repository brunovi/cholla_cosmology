#ifdef PARTICLES

#ifndef DENSITY_BOUNDARIES_H
#define DENSITY_BOUNDARIES_H


#include "particles_3D.h"

#ifdef MPI_CHOLLA
#include "../grid3D.h"
#endif


void Transfer_Particles_Density_Boundaries( Particles_3D &Parts );

#ifdef MPI_CHOLLA
void Transfer_Particles_Density_Boundaries_MPI( Grid3D &G, struct parameters P );
#endif





#endif
#endif
