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


void Advance_Particles( Grid3D &G );




#endif
#endif
