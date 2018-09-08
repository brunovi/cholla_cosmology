#ifdef PARTICLES

#ifndef DENSITY_CIC_CUDA_H
#define DENSITY_CIC_CUDA_H

#include<stdio.h>
#include<stdlib.h>
#include"math.h"
#include"../global.h"
#include"particles_3D.h"
#include <iostream>
// #include "../grid3D.h"




//
// void Clear_Density( Particles_3D &Parts );
//
// void Get_Indexes_CIC( Real xMin, Real yMin, Real zMin, Real dx, Real dy, Real dz, Real pos_x, Real pos_y, Real pos_z, int &indx_x, int &indx_y, int &indx_z );
//
// void Get_Particles_Density_CIC( Grid3D &G, struct parameters P, Real *time_pDens, Real *time_pDens_trans );
//
// void Copy_Particles_Density_to_Gravity( Grid3D &G );

// #ifdef PARTICLES_CUDA
void Get_Density_CIC_CUDA( Particles_3D &Parts );
// #endif

#endif
#endif
