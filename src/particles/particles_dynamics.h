#ifdef PARTICLES

#ifndef PARTICLES_DYNAMICS_H
#define PARTICLES_DYNAMICS_H

#include<stdio.h>
#include<stdlib.h>
#include"math.h"
#include"../global.h"
#include"../grid3D.h"
#include"particles_3D.h"
#include <iostream>
#include"gravity_CIC.h"

#ifdef PARTICLES_OMP
#include "particles_omp.h"
#endif

#ifdef COSMOLOGY
#include "../cosmology/cosmology.h"
#endif

// void Advance_Particles_LeapFrog( Particles_3D &Particles);

void Get_Particles_Acceleration( Grid3D &G, part_int_t p_start, part_int_t p_end, int g_start, int g_end );

void Advance_Particles_step1( Particles_3D &Particles, part_int_t p_start, part_int_t p_end );

void Advance_Particles_step2( Particles_3D &Particles, part_int_t p_start, part_int_t p_end );

#ifdef COSMOLOGY
void Advance_Particles_step1_cosmo( Particles_3D &Particles, Cosmology &Cosmo, part_int_t p_start, part_int_t p_end );
void Advance_Particles_step2_cosmo( Particles_3D &Particles, Cosmology &Cosmo, part_int_t p_start, part_int_t p_end );
#endif

Real Get_Particles_dt( Particles_3D &Particles );

Real Update_Particles( Grid3D &G, int step );






#endif
#endif
