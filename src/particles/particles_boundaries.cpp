#ifdef PARTICLES

#include "particles_boundaries.h"

void Tranfer_Particles_Boundaries( Particles_3D &Particles ){

  Real xMin, xMax, yMin, yMax, zMin, zMax, Lx, Ly, Lz;
  xMin = Particles.G.xMin;
  yMin = Particles.G.yMin;
  zMin = Particles.G.zMin;

  xMax = Particles.G.xMax;
  yMax = Particles.G.yMax;
  zMax = Particles.G.zMax;

  Lx = xMax - xMin;
  Ly = yMax - yMin;
  Lz = zMax - zMin;

  part_int_t pIndx;
  Real pPos_x, pPos_y, pPos_z;

  for ( pIndx=0; pIndx<Particles.n_local; pIndx++){
    pPos_x = Particles.pos_x[pIndx];
    if (pPos_x >= xMax) Particles.pos_x[pIndx] -= Lx;
    if (pPos_x <  xMin) Particles.pos_x[pIndx] += Lx;
  }

  for ( pIndx=0; pIndx<Particles.n_local; pIndx++){
    pPos_y = Particles.pos_y[pIndx];
    if (pPos_y >= yMax) Particles.pos_y[pIndx] -= Ly;
    if (pPos_y <  yMin) Particles.pos_y[pIndx] += Ly;
  }

  for ( pIndx=0; pIndx<Particles.n_local; pIndx++){
    pPos_z = Particles.pos_z[pIndx];
    if (pPos_z >= zMax) Particles.pos_z[pIndx] -= Lz;
    if (pPos_z <  zMin) Particles.pos_z[pIndx] += Lz;
  }

}

#endif
