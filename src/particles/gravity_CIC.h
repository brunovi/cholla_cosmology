#ifdef PARTICLES

#ifndef GRAVITY_CIC_H
#define GRAVITY_CIC_H

#include<stdio.h>
#include<stdlib.h>
#include"math.h"
#include"../global.h"
#include"../grid3D.h"
#include"particles_3D.h"
#include <iostream>
#include"density_CIC.h"


void Get_Gavity_Field( Grid3D &G );

void Get_Gravity_CIC( Grid3D &G );


#endif
#endif
