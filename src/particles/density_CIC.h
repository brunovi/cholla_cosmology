#ifdef PARTICLES

#ifndef DENSITY_CIC_H
#define DENSITY_CIC_H

#include<stdio.h>
#include<stdlib.h>
#include"math.h"
#include"../global.h"
#include"particles_3D.h"
#include <iostream>


void Clear_Density( Particles_3D &Parts );

void Get_Particles_Density_CIC( Particles_3D &Parts );

#endif
#endif
