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

// void Advance_Particles_LeapFrog( Particles_3D &Particles);

void Get_Particles_Acceleration( Grid3D &G );

void Advance_Particles_step1( Particles_3D &Particles );

void Advance_Particles_step2( Particles_3D &Particles );

float Get_Particles_dt( Particles_3D &Particles );

float Update_Particles( Grid3D &G, int step );






#endif
#endif
