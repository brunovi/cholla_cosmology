#ifdef COSMOLOGY

#ifndef COSMOLOGY_UNITS_H
#define COSMOLOGY_UNITS_H

#include<stdio.h>
#include <cmath>
#include "../io.h"
#include"../global.h"
#include"../grid3D.h"

void Change_DM_Frame_System( Grid3D &G, bool forward );
void Change_Cosmological_Frame_Sytem( Grid3D &G, bool forward );


#endif
#endif
