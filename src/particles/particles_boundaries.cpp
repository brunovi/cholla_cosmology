#ifdef PARTICLES

#include "particles_boundaries.h"
#include <unistd.h>
#include <algorithm>
#include <iostream>

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

#ifdef MPI_CHOLLA


int real_to_int( Real inVal ){
  int outVal = (int) inVal;
  if ( (inVal - outVal) > 0.1 ) outVal += 1;
  if ( abs(outVal - inVal) > 0.5 ) outVal -= 1;
  return outVal;
}


part_int_t real_to_part_int( Real inVal ){
  part_int_t outVal = (part_int_t) inVal;
  if ( (inVal - outVal) > 0.1 ) outVal += 1;
  if ( abs(outVal - inVal) > 0.5 ) outVal -= 1;
  return outVal;
}

Real Get_and_Remove_Real( part_int_t indx, real_vector_t &vec ){
  Real value = vec[indx];
  vec[indx] = vec.back();
  vec.pop_back();
  return value;
}

Real Get_and_Remove_partID( part_int_t indx, int_vector_t &vec ){
  Real value = (Real) vec[indx];
  vec[indx] = vec.back();
  vec.pop_back();
  return value;
}

void Remove_ID( part_int_t indx, int_vector_t &vec ){
  vec[indx] = vec.back();
  vec.pop_back();
}

void Remove_Real( part_int_t indx, real_vector_t &vec ){
  vec[indx] = vec.back();
  vec.pop_back();
}

void Particles_3D::Clear_Vectors_For_Transfers( void ){
  out_indxs_vec_x0.clear();
  out_indxs_vec_x1.clear();
  out_indxs_vec_y0.clear();
  out_indxs_vec_y1.clear();
  out_indxs_vec_z0.clear();
  out_indxs_vec_z1.clear();
}

void Particles_3D::Clear_Particles_For_Transfer( void ){
  n_transfer_x0 = 0;
  n_transfer_x1 = 0;
  n_transfer_y0 = 0;
  n_transfer_y1 = 0;
  n_transfer_z0 = 0;
  n_transfer_z1 = 0;

  n_send_x0 = 0;
  n_send_x1 = 0;
  n_send_y0 = 0;
  n_send_y1 = 0;
  n_send_z0 = 0;
  n_send_z1 = 0;

  n_recv_x0 = 0;
  n_recv_x1 = 0;
  n_recv_y0 = 0;
  n_recv_y1 = 0;
  n_recv_z0 = 0;
  n_recv_z1 = 0;

  n_in_buffer_x0 = 0;
  n_in_buffer_x1 = 0;
  n_in_buffer_y0 = 0;
  n_in_buffer_y1 = 0;
  n_in_buffer_z0 = 0;
  n_in_buffer_z1 = 0;
}


// void Particles_3D::Select_Particles_to_Transfer( int dir ){
//
//   out_indxs_vec_x0.clear();
//   out_indxs_vec_x1.clear();
//   out_indxs_vec_y0.clear();
//   out_indxs_vec_y1.clear();
//   out_indxs_vec_z0.clear();
//   out_indxs_vec_z1.clear();
//   //
//
//   part_int_t pIndx;
//   part_int_t pID;
//   for ( pIndx=0; pIndx<n_local; pIndx++ ){
//     if ( pos_x[pIndx] < G.xMin ){
//       out_indxs_vec_x0.push_back( pIndx );
//       continue;
//     }
//     if ( pos_x[pIndx] >= G.xMax ){
//       out_indxs_vec_x1.push_back( pIndx );
//       continue;
//     }
//     if ( pos_y[pIndx] < G.yMin ){
//       out_indxs_vec_y0.push_back( pIndx );
//       continue;
//     }
//     if ( pos_y[pIndx] >= G.yMax ){
//       out_indxs_vec_y1.push_back( pIndx );
//       continue;
//     }
//     if ( pos_z[pIndx] < G.zMin ){
//       out_indxs_vec_z0.push_back( pIndx );
//       continue;
//     }
//     if ( pos_z[pIndx] >= G.zMax ){
//       out_indxs_vec_z1.push_back( pIndx );
//       continue;
//     }
//   }
//   std::sort(out_indxs_vec_x0.begin(), out_indxs_vec_x0.end());
//   std::sort(out_indxs_vec_x1.begin(), out_indxs_vec_x1.end());
//   std::sort(out_indxs_vec_y0.begin(), out_indxs_vec_y0.end());
//   std::sort(out_indxs_vec_y1.begin(), out_indxs_vec_y1.end());
//   std::sort(out_indxs_vec_z0.begin(), out_indxs_vec_z0.end());
//   std::sort(out_indxs_vec_z1.begin(), out_indxs_vec_z1.end());
// }




void Particles_3D::Select_Particles_to_Transfer_All(  ){
  //
  int i;
  bool already_transfered, pID_in_transfers;
  part_int_t pIndx;
  // part_int_t pID;
  // if ( dir == 0 ){
  out_indxs_vec_x0.clear();
  out_indxs_vec_x1.clear();
  out_indxs_vec_y0.clear();
  out_indxs_vec_y1.clear();
  out_indxs_vec_z0.clear();
  out_indxs_vec_z1.clear();

  for ( pIndx=0; pIndx<n_local; pIndx++ ){
    // pID = partIDs[pIndx];
    if ( pos_x[pIndx] < G.xMin ){
      out_indxs_vec_x0.push_back( pIndx );
      continue;
    }
    if ( pos_x[pIndx] >= G.xMax ){
      out_indxs_vec_x1.push_back( pIndx );
      continue;
    }
    if ( pos_y[pIndx] < G.yMin ){
      out_indxs_vec_y0.push_back( pIndx );
      continue;
    }
    if ( pos_y[pIndx] >= G.yMax ){
      out_indxs_vec_y1.push_back( pIndx );
      continue;
    }
    if ( pos_z[pIndx] < G.zMin ){
      out_indxs_vec_z0.push_back( pIndx );
      continue;
    }
    if ( pos_z[pIndx] >= G.zMax ){
      out_indxs_vec_z1.push_back( pIndx );
      continue;
    }
  }

  std::sort(out_indxs_vec_x0.begin(), out_indxs_vec_x0.end());
  std::sort(out_indxs_vec_x1.begin(), out_indxs_vec_x1.end());
  std::sort(out_indxs_vec_y0.begin(), out_indxs_vec_y0.end());
  std::sort(out_indxs_vec_y1.begin(), out_indxs_vec_y1.end());
  std::sort(out_indxs_vec_z0.begin(), out_indxs_vec_z0.end());
  std::sort(out_indxs_vec_z1.begin(), out_indxs_vec_z1.end());
  n_send_x0 += out_indxs_vec_x0.size();
  n_send_x1 += out_indxs_vec_x1.size();
  n_send_y0 += out_indxs_vec_y0.size();
  n_send_y1 += out_indxs_vec_y1.size();
  n_send_z0 += out_indxs_vec_z0.size();
  n_send_z1 += out_indxs_vec_z1.size();
}


void Particles_3D::Select_Particles_to_Transfer( int dir ){

  int i;
  bool already_transfered, pID_in_transfers;
  part_int_t pIndx;
  // part_int_t pID;
  if ( dir == 0 ){
    out_indxs_vec_x0.clear();
    out_indxs_vec_x1.clear();
    for ( pIndx=0; pIndx<n_local; pIndx++ ){
      // pID = partIDs[pIndx];
      if ( pos_x[pIndx] < G.xMin ){
        out_indxs_vec_x0.push_back( pIndx );
        continue;
      }
      if ( pos_x[pIndx] >= G.xMax ){
        out_indxs_vec_x1.push_back( pIndx );
        continue;
      }
    }
    std::sort(out_indxs_vec_x0.begin(), out_indxs_vec_x0.end());
    std::sort(out_indxs_vec_x1.begin(), out_indxs_vec_x1.end());
    n_send_x0 += out_indxs_vec_x0.size();
    n_send_x1 += out_indxs_vec_x1.size();
    return;
  }
  if ( dir == 1 ){
    out_indxs_vec_y0.clear();
    out_indxs_vec_y1.clear();
    for ( pIndx=0; pIndx<n_local; pIndx++ ){
      // pID = partIDs[pIndx];
      already_transfered = false;
      if ( pos_x[pIndx] < G.xMin || pos_x[pIndx] >= G.xMax ) already_transfered = true;
      if ( already_transfered ){
        continue;
      }
      if ( pos_y[pIndx] < G.yMin ){
        out_indxs_vec_y0.push_back( pIndx );
        continue;
      }
      if ( pos_y[pIndx] >= G.yMax ){
        out_indxs_vec_y1.push_back( pIndx );
        continue;
      }
    }
    std::sort(out_indxs_vec_y0.begin(), out_indxs_vec_y0.end());
    std::sort(out_indxs_vec_y1.begin(), out_indxs_vec_y1.end());
    n_send_y0 += out_indxs_vec_y0.size();
    n_send_y1 += out_indxs_vec_y1.size();
  }
  if ( dir == 2 ){
    out_indxs_vec_z0.clear();
    out_indxs_vec_z1.clear();
    for ( pIndx=0; pIndx<n_local; pIndx++ ){
      // pID = partIDs[pIndx];
      already_transfered = false;
      if ( pos_x[pIndx] < G.xMin || pos_x[pIndx] >= G.xMax || pos_y[pIndx] < G.yMin || pos_y[pIndx] >= G.yMax ) already_transfered = true;
      if ( already_transfered ){
        continue;
      }
      if ( pos_z[pIndx] < G.zMin ){
        out_indxs_vec_z0.push_back( pIndx );
        continue;
      }
      if ( pos_z[pIndx] >= G.zMax ){
        out_indxs_vec_z1.push_back( pIndx );
        continue;
      }
    }
    std::sort(out_indxs_vec_z0.begin(), out_indxs_vec_z0.end());
    std::sort(out_indxs_vec_z1.begin(), out_indxs_vec_z1.end());
    n_send_z0 += out_indxs_vec_z0.size();
    n_send_z1 += out_indxs_vec_z1.size();
  }
}



void Particles_3D::Load_Particles_to_Buffer_new( int direction, int side, Real *send_buffer, int buffer_length  ){

  part_int_t n_out;
  part_int_t n_send;
  int_vector_t *out_indxs_vec;
  part_int_t *n_in_buffer;

  if ( direction == 0 ){
    if ( side == 0 ){
      out_indxs_vec = &out_indxs_vec_x0;
      n_send = n_send_x0;
      n_in_buffer = &n_in_buffer_x0;
    }
    if ( side == 1 ){
      out_indxs_vec = &out_indxs_vec_x1;
      n_send = n_send_x1;
      n_in_buffer = &n_in_buffer_x1;
    }
  }
  if ( direction == 1 ){
    if ( side == 0 ){
      out_indxs_vec = &out_indxs_vec_y0;
      n_send = n_send_y0;
      n_in_buffer = &n_in_buffer_y0;
    }
    if ( side == 1 ){
      out_indxs_vec = &out_indxs_vec_y1;
      n_send = n_send_y1;
      n_in_buffer = &n_in_buffer_y1;
    }
  }
  if ( direction == 2 ){
    if ( side == 0 ){
      out_indxs_vec = &out_indxs_vec_z0;
      n_send = n_send_z0;
      n_in_buffer = &n_in_buffer_z0;
    }
    if ( side == 1 ){
      out_indxs_vec = &out_indxs_vec_z1;
      n_send = n_send_z1;
      n_in_buffer = &n_in_buffer_z1;
    }
  }

  part_int_t offset, offset_extra;

  n_out = out_indxs_vec->size();
  offset = *n_in_buffer*N_DATA_PER_PARTICLE_TRANSFER;

  part_int_t indx, pIndx;
  for ( indx=0; indx<n_out; indx++ ){
    pIndx = out_indxs_vec->back();
    send_buffer[ offset + 0 ] = Get_and_Remove_Real( pIndx, pos_x );
    send_buffer[ offset + 1 ] = Get_and_Remove_Real( pIndx, pos_y );
    send_buffer[ offset + 2 ] = Get_and_Remove_Real( pIndx, pos_z );
    send_buffer[ offset + 3 ] = Get_and_Remove_Real( pIndx, vel_x );
    send_buffer[ offset + 4 ] = Get_and_Remove_Real( pIndx, vel_y );
    send_buffer[ offset + 5 ] = Get_and_Remove_Real( pIndx, vel_z );

    offset_extra = offset + 5;
    #ifndef SINGLE_PARTICLE_MASS
    offset_extra += 1;
    send_buffer[ offset_extra ] = Get_and_Remove_Real( pIndx, mass );
    #endif
    #ifdef PARTICLE_IDS
    offset_extra += 1;
    send_buffer[ offset_extra ] = (Real) Get_and_Remove_partID( pIndx, partIDs );
    #endif

    Remove_Real( pIndx, grav_x );
    Remove_Real( pIndx, grav_y );
    Remove_Real( pIndx, grav_z );
    // *n_send += 1;
    n_local -= 1;
    *n_in_buffer += 1;
    out_indxs_vec->pop_back();
    offset += N_DATA_PER_PARTICLE_TRANSFER;
    if ( offset > buffer_length ) std::cout << "ERROR: Buffer length exceeded on particles transfer" << std::endl;
  }

  if (out_indxs_vec->size() > 0 ) std::cout << "ERROR: Particle output vector not empty after transfer " << std::endl;
}

//
// void Particles_3D::Load_Particles_to_Buffer( int direction, int side, int buffer_start, Real *send_buffer, int MAX_PARTICLES_IN_BUFFER, bool secondary  ){
//
//   part_int_t n_out;
//   part_int_t *n_send;
//   int_vector_t *out_indxs_vec;
//
//   if ( direction == 0 ){
//     if ( side == 0 ){
//       out_indxs_vec = &out_indxs_vec_x0;
//       n_send = &n_send_x0;
//     }
//     if ( side == 1 ){
//       out_indxs_vec = &out_indxs_vec_x1;
//       n_send = &n_send_x1;
//     }
//   }
//   if ( direction == 1 ){
//     if ( side == 0 ){
//       out_indxs_vec = &out_indxs_vec_y0;
//       n_send = &n_send_y0;
//     }
//     if ( side == 1 ){
//       out_indxs_vec = &out_indxs_vec_y1;
//       n_send = &n_send_y1;
//     }
//   }
//   if ( direction == 2 ){
//     if ( side == 0 ){
//       out_indxs_vec = &out_indxs_vec_z0;
//       n_send = &n_send_z0;
//     }
//     if ( side == 1 ){
//       out_indxs_vec = &out_indxs_vec_z1;
//       n_send = &n_send_z1;
//     }
//   }
//
//   part_int_t offset, n_in_buffer, offset_extra;
//   part_int_t indx, pIndx, indx_start;
//
//   // if (!secondary){
//   n_out = out_indxs_vec->size();
//   n_in_buffer = real_to_part_int( send_buffer[buffer_start] ) ;
//   if ( (n_out + n_in_buffer) >= MAX_PARTICLES_IN_BUFFER ) n_out = MAX_PARTICLES_IN_BUFFER - n_in_buffer;
//   // *n_transfer_1 = n_out;
//   // *n_transfer_2 = out_indxs_vec->size() - n_out;
//   // indx_start = 0;
// // }
// // else{
//   // indx_start = *n_transfer_1;
//   // n_out = *n_transfer_2;
//
//   // }
//
//   offset = buffer_start + N_HEADER_PARTICLES_TRANSFER + n_in_buffer*N_DATA_PER_PARTICLE_TRANSFER;
//   for ( indx=0; indx<n_out; indx++ ){
//     pIndx = out_indxs_vec->back();
//     send_buffer[ offset + 0 ] = Get_and_Remove_Real( pIndx, pos_x );
//     send_buffer[ offset + 1 ] = Get_and_Remove_Real( pIndx, pos_y );
//     send_buffer[ offset + 2 ] = Get_and_Remove_Real( pIndx, pos_z );
//     send_buffer[ offset + 3 ] = Get_and_Remove_Real( pIndx, vel_x );
//     send_buffer[ offset + 4 ] = Get_and_Remove_Real( pIndx, vel_y );
//     send_buffer[ offset + 5 ] = Get_and_Remove_Real( pIndx, vel_z );
//
//     offset_extra = offset + 5;
//     #ifndef SINGLE_PARTICLE_MASS
//     offset_extra += 1;
//     send_buffer[ offset_extra ] = Get_and_Remove_Real( pIndx, mass );
//     #endif
//     #ifdef PARTICLE_IDS
//     offset_extra += 1;
//     send_buffer[ offset_extra ] = (Real) Get_and_Remove_partID( pIndx, partIDs );
//     #endif
//
//     Remove_Real( pIndx, grav_x );
//     Remove_Real( pIndx, grav_y );
//     Remove_Real( pIndx, grav_z );
//
//     send_buffer[buffer_start] += 1;
//     // *n_send += 1;
//     n_local -= 1;
//     out_indxs_vec->pop_back();
//     offset += N_DATA_PER_PARTICLE_TRANSFER;
//   }
//
//   if (secondary && (out_indxs_vec->size() > 0 ) ) std::cout << "ERROR: Particle output vector not empty after secondary transfer " << std::endl;
//   send_buffer[buffer_start+1] += out_indxs_vec->size();
//
//
// }


void Particles_3D::Add_Particle_To_Buffer_new( Real *buffer, part_int_t n_in_buffer, int buffer_length, Real pId, Real pMass,
                            Real pPos_x, Real pPos_y, Real pPos_z, Real pVel_x, Real pVel_y, Real pVel_z){

  int offset, offset_extra;
  offset = n_in_buffer * N_DATA_PER_PARTICLE_TRANSFER;

  if (offset > buffer_length ) if ( offset > buffer_length ) std::cout << "ERROR: Buffer length exceeded on particles transfer" << std::endl;
  buffer[offset + 0] = pPos_x;
  buffer[offset + 1] = pPos_y;
  buffer[offset + 2] = pPos_z;
  buffer[offset + 3] = pVel_x;
  buffer[offset + 4] = pVel_y;
  buffer[offset + 5] = pVel_z;

  offset_extra = offset + 5;
  #ifndef SINGLE_PARTICLE_MASS
  offset_extra += 1;
  buffer[ offset_extra ] = pMass;
  #endif
  #ifdef PARTICLE_IDS
  offset_extra += 1;
  buffer[offset_extra] = pId;
  #endif


}


// void Particles_3D::Add_Particle_To_Buffer( Real *buffer, int buffer_start, int max_particles, Real *buffer_secondary, int buffer_start_secondary, Real pId, Real pMass,
//                             Real pPos_x, Real pPos_y, Real pPos_z,
//                             Real pVel_x, Real pVel_y, Real pVel_z, bool secondary ){
//   int n_in_buffer, offset_out, offset_extra;
//   if ( !secondary ){
//     n_in_buffer = real_to_part_int( buffer[buffer_start] );
//     if ( n_in_buffer < max_particles ){
//       offset_out = buffer_start + N_HEADER_PARTICLES_TRANSFER + n_in_buffer*N_DATA_PER_PARTICLE_TRANSFER;
//       buffer[offset_out + 0] = pPos_x;
//       buffer[offset_out + 1] = pPos_y;
//       buffer[offset_out + 2] = pPos_z;
//       buffer[offset_out + 3] = pVel_x;
//       buffer[offset_out + 4] = pVel_y;
//       buffer[offset_out + 5] = pVel_z;
//
//       offset_extra = offset_out + 5;
//       #ifndef SINGLE_PARTICLE_MASS
//       offset_extra += 1;
//       buffer[ offset_extra ] = pMass;
//       #endif
//       #ifdef PARTICLE_IDS
//       offset_extra += 1;
//       buffer[offset_extra] = pId;
//       #endif
//
//       buffer[buffer_start] += 1;
//     }
//     else{
//       buffer[buffer_start+1] += 1;
//       n_in_buffer = real_to_part_int( buffer_secondary[buffer_start_secondary] );
//       offset_out = buffer_start_secondary + N_HEADER_PARTICLES_TRANSFER + n_in_buffer*N_DATA_PER_PARTICLE_TRANSFER;
//       buffer_secondary[offset_out + 0] = pPos_x;
//       buffer_secondary[offset_out + 1] = pPos_y;
//       buffer_secondary[offset_out + 2] = pPos_z;
//       buffer_secondary[offset_out + 3] = pVel_x;
//       buffer_secondary[offset_out + 4] = pVel_y;
//       buffer_secondary[offset_out + 5] = pVel_z;
//
//       offset_extra = offset_out + 5;
//       #ifndef SINGLE_PARTICLE_MASS
//       offset_extra += 1;
//       buffer_secondary[ offset_extra ] = pMass;
//       #endif
//       #ifdef PARTICLE_IDS
//       offset_extra += 1;
//       buffer_secondary[offset_extra] = pId;
//       #endif
//
//       buffer_secondary[buffer_start_secondary] += 1;
//     }
//   }
//   else{
//     n_in_buffer = real_to_part_int( buffer[buffer_start] );
//     offset_out = buffer_start + N_HEADER_PARTICLES_TRANSFER + n_in_buffer*N_DATA_PER_PARTICLE_TRANSFER;
//
//     buffer[offset_out + 0] = pPos_x;
//     buffer[offset_out + 1] = pPos_y;
//     buffer[offset_out + 2] = pPos_z;
//     buffer[offset_out + 3] = pVel_x;
//     buffer[offset_out + 4] = pVel_y;
//     buffer[offset_out + 5] = pVel_z;
//
//     offset_extra = offset_out + 5;
//     #ifndef SINGLE_PARTICLE_MASS
//     offset_extra += 1;
//     buffer[offset_extra ] = pMass;
//     #endif
//     #ifdef PARTICLE_IDS
//     offset_extra += 1;
//     buffer[offset_extra] = pId;
//     #endif
//
//     buffer[buffer_start] += 1;
//   }
// }

void Particles_3D::Add_Particle_To_Vectors( Real pId, Real pMass,
                            Real pPos_x, Real pPos_y, Real pPos_z,
                            Real pVel_x, Real pVel_y, Real pVel_z ){
  pos_x.push_back( pPos_x );
  pos_y.push_back( pPos_y );
  pos_z.push_back( pPos_z );
  vel_x.push_back( pVel_x );
  vel_y.push_back( pVel_y );
  vel_z.push_back( pVel_z );
  #ifndef SINGLE_PARTICLE_MASS
  mass.push_back( pMass );
  #endif
  #ifdef PARTICLE_IDS
  partIDs.push_back( real_to_part_int(pId) );
  #endif
  grav_x.push_back(0);
  grav_y.push_back(0);
  grav_z.push_back(0);

  n_local += 1;
}


void Particles_3D::Unload_Particles_from_Buffer_new( int direction, int side, Real *recv_buffer, part_int_t n_recv,
      Real *send_buffer_y0, Real *send_buffer_y1, Real *send_buffer_z0, Real *send_buffer_z1, int buffer_length_y0, int buffer_length_y1, int buffer_length_z0, int buffer_length_z1){

  int offset_buff, offset_extra;
  part_int_t pId;
  Real pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z;

  offset_buff = 0;
  part_int_t indx;
  for ( indx=0; indx<n_recv; indx++ ){
    pPos_x = recv_buffer[ offset_buff + 0 ];
    pPos_y = recv_buffer[ offset_buff + 1 ];
    pPos_z = recv_buffer[ offset_buff + 2 ];
    pVel_x = recv_buffer[ offset_buff + 3 ];
    pVel_y = recv_buffer[ offset_buff + 4 ];
    pVel_z = recv_buffer[ offset_buff + 5 ];

    offset_extra = offset_buff + 5;
    #if SINGLE_PARTICLE_MASS
    pMass = particle_mass;
    #else
    offset_extra += 1;
    pMass  = recv_buffer[ offset_extra ];
    #endif
    #ifdef PARTICLE_IDS
    offset_extra += 1;
    pId    = recv_buffer[ offset_extra ];
    #else
    pId = 0;
    #endif

    offset_buff += N_DATA_PER_PARTICLE_TRANSFER;
    if ( pPos_x <  G.domainMin_x ) pPos_x += ( G.domainMax_x - G.domainMin_x );
    if ( pPos_x >= G.domainMax_x ) pPos_x -= ( G.domainMax_x - G.domainMin_x );
    if ( ( pPos_x < G.xMin ) || ( pPos_x >= G.xMax )  ){
      #ifdef PARTICLE_IDS
      std::cout << "ERROR Particle Transfer out of X domain    pID: " << pId << std::endl;
      #else
      std::cout << "ERROR Particle Transfer out of X domain" << std::endl;
      #endif
      std::cout << " posX: " << pPos_x << " velX: " << pVel_x << std::endl;
      std::cout << " posY: " << pPos_y << " velY: " << pVel_y << std::endl;
      std::cout << " posZ: " << pPos_z << " velZ: " << pVel_z << std::endl;
      std::cout << " Domain X: " << G.xMin << "  " << G.xMax << std::endl;
      std::cout << " Domain Y: " << G.yMin << "  " << G.yMax << std::endl;
      std::cout << " Domain Z: " << G.zMin << "  " << G.zMax << std::endl;
      continue;
    }

    if (  direction == 1 ){
      if ( pPos_y <  G.domainMin_y ) pPos_y += ( G.domainMax_y - G.domainMin_y );
      if ( pPos_y >= G.domainMax_y ) pPos_y -= ( G.domainMax_y - G.domainMin_y );
    }
    if ( (direction==1 || direction==2) && (( pPos_y < G.yMin ) || ( pPos_y >= G.yMax ))  ){
      #ifdef PARTICLE_IDS
      std::cout << "ERROR Particle Transfer out of Y domain    pID: " << pId << std::endl;
      #else
      std::cout << "ERROR Particle Transfer out of Y domain" << std::endl;
      #endif
      std::cout << " posX: " << pPos_x << " velX: " << pVel_x << std::endl;
      std::cout << " posY: " << pPos_y << " velY: " << pVel_y << std::endl;
      std::cout << " posZ: " << pPos_z << " velZ: " << pVel_z << std::endl;
      std::cout << " Domain X: " << G.xMin << "  " << G.xMax << std::endl;
      std::cout << " Domain Y: " << G.yMin << "  " << G.yMax << std::endl;
      std::cout << " Domain Z: " << G.zMin << "  " << G.zMax << std::endl;
      continue;
    }

    if (direction  == 0 ){
      if ( pPos_y < G.yMin ){
        // std::cout << "Added Y0" << std::endl;
        Add_Particle_To_Buffer_new( send_buffer_y0, n_in_buffer_y0, buffer_length_y0, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
        n_send_y0 += 1;
        n_in_buffer_y0 += 1;
        continue;
      }
      if ( pPos_y >= G.yMax ){
        // std::cout << "Added Y1" << std::endl;
        Add_Particle_To_Buffer_new( send_buffer_y1, n_in_buffer_y1, buffer_length_y1, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
        n_send_y1 += 1;
        n_in_buffer_y1 += 1;
        continue;
      }
    }

    if (  direction == 2 ){
      if ( pPos_z <  G.domainMin_z ) pPos_z += ( G.domainMax_z - G.domainMin_z );
      if ( pPos_z >= G.domainMax_z ) pPos_z -= ( G.domainMax_z - G.domainMin_z );
    }
    if ( (direction==2) && (( pPos_z < G.zMin ) || ( pPos_z >= G.zMax ))  ){
      #ifdef PARTICLE_IDS
      std::cout << "ERROR Particle Transfer out of Z domain    pID: " << pId << std::endl;
      #else
      std::cout << "ERROR Particle Transfer out of Z domain" << std::endl;
      #endif
      std::cout << " posX: " << pPos_x << " velX: " << pVel_x << std::endl;
      std::cout << " posY: " << pPos_y << " velY: " << pVel_y << std::endl;
      std::cout << " posZ: " << pPos_z << " velZ: " << pVel_z << std::endl;
      std::cout << " Domain X: " << G.xMin << "  " << G.xMax << std::endl;
      std::cout << " Domain Y: " << G.yMin << "  " << G.yMax << std::endl;
      std::cout << " Domain Z: " << G.zMin << "  " << G.zMax << std::endl;
      continue;
    }
  //
    if (direction  !=2 ){
      if ( pPos_z < G.zMin ){
        // std::cout << "Added Z0" << std::endl;
        Add_Particle_To_Buffer_new( send_buffer_z0, n_in_buffer_z0, buffer_length_z0, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
        n_send_z0 += 1;
        n_in_buffer_z0 += 1;
        continue;
      }
      if ( pPos_z >= G.zMax ){
        // std::cout << "Added Z1" << std::endl;
        Add_Particle_To_Buffer_new( send_buffer_z1, n_in_buffer_z1, buffer_length_z1, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
        n_send_z1 += 1;
        n_in_buffer_z1 += 1;
        continue;
      }
    }
    Add_Particle_To_Vectors( pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
  }
}


//
// void Particles_3D::Unload_Particles_from_Buffer( int direction, int side, int buffer_start, Real *recv_buffer,
//       Real *send_buffer_y0, Real *send_buffer_y1, Real *send_buffer_z0, Real *send_buffer_z1, int buffer_start_y, int buffer_start_z,
//       Real *send_buffer_y0_second, Real *send_buffer_y1_second, Real *send_buffer_z0_second, Real *send_buffer_z1_second, int max_particles, bool secondary  ){
//
//   int n_recv;
//   int offset_buff, indx, offset_out, n_in_out, offset_extra;
//   part_int_t pId;
//   Real pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z;
//
//   n_recv = real_to_part_int( recv_buffer[buffer_start] );
//   offset_buff = buffer_start + N_HEADER_PARTICLES_TRANSFER;
//   for ( indx = 0; indx<n_recv; indx++ ){
//     pPos_x = recv_buffer[ offset_buff + 0 ];
//     pPos_y = recv_buffer[ offset_buff + 1 ];
//     pPos_z = recv_buffer[ offset_buff + 2 ];
//     pVel_x = recv_buffer[ offset_buff + 3 ];
//     pVel_y = recv_buffer[ offset_buff + 4 ];
//     pVel_z = recv_buffer[ offset_buff + 5 ];
//
//     offset_extra = offset_buff + 5;
//     #if SINGLE_PARTICLE_MASS
//     pMass = particle_mass;
//     #else
//     offset_extra += 1;
//     pMass  = recv_buffer[ offset_extra ];
//     #endif
//     #ifdef PARTICLE_IDS
//     offset_extra += 1;
//     pId    = recv_buffer[ offset_extra ];
//     #else
//     pId = 0;
//     #endif
//
//
//     offset_buff += N_DATA_PER_PARTICLE_TRANSFER;
//     if ( pPos_x <  G.domainMin_x ) pPos_x += ( G.domainMax_x - G.domainMin_x );
//     if ( pPos_x >= G.domainMax_x ) pPos_x -= ( G.domainMax_x - G.domainMin_x );
//     if ( ( pPos_x < G.xMin ) || ( pPos_x >= G.xMax )  ){
//       #ifdef PARTICLE_IDS
//       std::cout << "ERROR Particle Transfer out of X domain    pID: " << pId << std::endl;
//       #else
//       std::cout << "ERROR Particle Transfer out of X domain" << std::endl;
//       #endif
//       std::cout << " posX: " << pPos_x << " velX: " << pVel_x << std::endl;
//       std::cout << " posY: " << pPos_y << " velY: " << pVel_y << std::endl;
//       std::cout << " posZ: " << pPos_z << " velZ: " << pVel_z << std::endl;
//       std::cout << " Domain X: " << G.xMin << "  " << G.xMax << std::endl;
//       std::cout << " Domain Y: " << G.yMin << "  " << G.yMax << std::endl;
//       std::cout << " Domain Z: " << G.zMin << "  " << G.zMax << std::endl;
//
//       continue;
//     }
//
//     if (  direction == 1 ){
//       if ( pPos_y <  G.domainMin_y ) pPos_y += ( G.domainMax_y - G.domainMin_y );
//       if ( pPos_y >= G.domainMax_y ) pPos_y -= ( G.domainMax_y - G.domainMin_y );
//     }
//     if ( (direction==1 || direction==2) && (( pPos_y < G.yMin ) || ( pPos_y >= G.yMax ))  ){
//       #ifdef PARTICLE_IDS
//       std::cout << "ERROR Particle Transfer out of Y domain    pID: " << pId << std::endl;
//       #else
//       std::cout << "ERROR Particle Transfer out of Y domain" << std::endl;
//       #endif
//       std::cout << " posX: " << pPos_x << " velX: " << pVel_x << std::endl;
//       std::cout << " posY: " << pPos_y << " velY: " << pVel_y << std::endl;
//       std::cout << " posZ: " << pPos_z << " velZ: " << pVel_z << std::endl;
//       std::cout << " Domain X: " << G.xMin << "  " << G.xMax << std::endl;
//       std::cout << " Domain Y: " << G.yMin << "  " << G.yMax << std::endl;
//       std::cout << " Domain Z: " << G.zMin << "  " << G.zMax << std::endl;
//       continue;
//     }
//
//     if (direction  == 0 ){
//       if ( pPos_y < G.yMin ){
//         // n_send_y0 += 1;
//         Add_Particle_To_Buffer( send_buffer_y0, buffer_start_y, max_particles, send_buffer_y0_second, 0, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z, secondary);
//         continue;
//       }
//       if ( pPos_y >= G.yMax ){
//         // n_send_y1 += 1;
//         Add_Particle_To_Buffer( send_buffer_y1, buffer_start_y, max_particles, send_buffer_y1_second, 0, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z, secondary);
//         continue;
//       }
//     }
//
//     if (  direction == 2 ){
//       if ( pPos_z <  G.domainMin_z ) pPos_z += ( G.domainMax_z - G.domainMin_z );
//       if ( pPos_z >= G.domainMax_z ) pPos_z -= ( G.domainMax_z - G.domainMin_z );
//     }
//     if ( (direction==2) && (( pPos_z < G.zMin ) || ( pPos_z >= G.zMax ))  ){
//       #ifdef PARTICLE_IDS
//       std::cout << "ERROR Particle Transfer out of Z domain    pID: " << pId << std::endl;
//       #else
//       std::cout << "ERROR Particle Transfer out of Z domain" << std::endl;
//       #endif
//       std::cout << " posX: " << pPos_x << " velX: " << pVel_x << std::endl;
//       std::cout << " posY: " << pPos_y << " velY: " << pVel_y << std::endl;
//       std::cout << " posZ: " << pPos_z << " velZ: " << pVel_z << std::endl;
//       std::cout << " Domain X: " << G.xMin << "  " << G.xMax << std::endl;
//       std::cout << " Domain Y: " << G.yMin << "  " << G.yMax << std::endl;
//       std::cout << " Domain Z: " << G.zMin << "  " << G.zMax << std::endl;
//       continue;
//     }
//   //
//     if (direction  !=2 ){
//       if ( pPos_z < G.zMin ){
//         n_send_z0 += 1;
//         Add_Particle_To_Buffer( send_buffer_z0, buffer_start_z, max_particles, send_buffer_z0_second, 0, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z, secondary);
//         continue;
//       }
//       if ( pPos_z >= G.zMax ){
//         n_send_z1 += 1;
//         Add_Particle_To_Buffer( send_buffer_z1, buffer_start_z, max_particles, send_buffer_z1_second, 0, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z, secondary);
//         continue;
//       }
//     }
//     Add_Particle_To_Vectors( pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
//   }
// }

void Particles_3D::Remove_Transfered_Particles( void ){
  int n_in_out_vectors, n_in_vectors;

  n_in_vectors =  pos_x.size() + pos_y.size() + pos_z.size() + vel_x.size() + vel_y.size() + vel_z.size() ;
  #ifndef SINGLE_PARTICLE_MASS
  n_in_vectors += mass.size();
  #endif
  #ifdef PARTICLE_IDS
  n_in_vectors += partIDs.size();
  #endif



  if ( n_in_vectors != n_local * N_DATA_PER_PARTICLE_TRANSFER ){
    std::cout << "ERROR PARTICLES TRANSFER: DATA IN VECTORS DIFFERENT FROM N_LOCAL###########" << std::endl;
    exit(-1);
  }

  n_in_out_vectors = out_indxs_vec_x0.size() + out_indxs_vec_x1.size() + out_indxs_vec_y0.size() + out_indxs_vec_y1.size() + out_indxs_vec_z0.size() + out_indxs_vec_z1.size();
  if ( n_in_out_vectors != 0 ){
    std::cout << "#################ERROR PARTICLES TRANSFER: OUPTUT VECTORS NOT EMPTY, N_IN_VECTORS: " << n_in_out_vectors << std::endl;
    part_int_t pId;
    if ( out_indxs_vec_x0.size()>0){
      std::cout << " In x0" << std::endl;
      pId = out_indxs_vec_x0[0];
    }
    if ( out_indxs_vec_x1.size()>0){
      std::cout << " In x1" << std::endl;
      pId = out_indxs_vec_x1[0];
    }
    if ( out_indxs_vec_y0.size()>0){
      std::cout << " In y0" << std::endl;
      pId = out_indxs_vec_y0[0];
    }
    if ( out_indxs_vec_y1.size()>0){
      std::cout << " In y1" << std::endl;
      pId = out_indxs_vec_y1[0];
    }
    if ( out_indxs_vec_z0.size()>0){
      std::cout << " In z0" << std::endl;
      pId = out_indxs_vec_z0[0];
    }
    if ( out_indxs_vec_z1.size()>0){
      std::cout << " In z1" << std::endl;
      pId = out_indxs_vec_z1[0];
    }
    std::cout  << "pos_x: " << pos_x[pId] << " x: " << G.xMin << "  " << G.xMax << std::endl;
    std::cout  << "pos_y: " << pos_y[pId] << " y: " << G.yMin << "  " << G.yMax << std::endl;
    std::cout  << "pos_z: " << pos_z[pId] << " z: " << G.zMin << "  " << G.zMax << std::endl;
    exit(-1);
  }

}


#endif
#endif
