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

void Particles_3D::Select_Particles_to_Transfer( void ){

  Clear_Vectors_For_Transfers();

  part_int_t pIndx;
  for ( pIndx=0; pIndx<n_local; pIndx++ ){
    // if (procID == 0) std::cout << pIndx << std::endl;
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

}

void Particles_3D::Load_Particles_to_Buffer( int direction, int side, int buffer_start, Real *send_buffer,int MAX_PARTICLES_IN_BUFFER  ){

  int n_out;
  int_vector_t *out_indxs_vec;

  if ( direction == 0 ){
    if ( side == 0 ){
      out_indxs_vec = &out_indxs_vec_x0;
    }
    if ( side == 1 ){
      out_indxs_vec = &out_indxs_vec_x1;
    }
  }
  if ( direction == 1 ){
    if ( side == 0 ){
      out_indxs_vec = &out_indxs_vec_y0;
    }
    if ( side == 1 ){
      out_indxs_vec = &out_indxs_vec_y1;
    }
  }
  if ( direction == 2 ){
    if ( side == 0 ){
      out_indxs_vec = &out_indxs_vec_z0;
    }
    if ( side == 1 ){
      out_indxs_vec = &out_indxs_vec_z1;
    }
  }

  n_out = out_indxs_vec->size();

  int offset, n_in_buffer;
  part_int_t indx, pIndx;

  n_in_buffer = real_to_int( send_buffer[buffer_start] ) ;

  offset = buffer_start + N_HEADER_PARTICLES_TRANSFER + n_in_buffer*N_DATA_PER_PARTICLE_TRANSFER;
  for ( indx=0; indx<n_out; indx++ ){
    pIndx = out_indxs_vec->back();
    send_buffer[ offset + 0 ] = (Real) Get_and_Remove_partID( pIndx, partIDs );
    send_buffer[ offset + 1 ] = Get_and_Remove_Real( pIndx, mass );
    send_buffer[ offset + 2 ] = Get_and_Remove_Real( pIndx, pos_x );
    send_buffer[ offset + 3 ] = Get_and_Remove_Real( pIndx, pos_y );
    send_buffer[ offset + 4 ] = Get_and_Remove_Real( pIndx, pos_z );
    send_buffer[ offset + 5 ] = Get_and_Remove_Real( pIndx, vel_x );
    send_buffer[ offset + 6 ] = Get_and_Remove_Real( pIndx, vel_y );
    send_buffer[ offset + 7 ] = Get_and_Remove_Real( pIndx, vel_z );
    send_buffer[buffer_start] += 1;
    n_local -= 1;
    out_indxs_vec->pop_back();
    offset += N_DATA_PER_PARTICLE_TRANSFER;
  }
  send_buffer[buffer_start+1] = 0;
  // for ( indx = n_out-1; indx>=0; indx-- ){
  //   pIndx = (*out_indxs_vec)[indx];
  //   // pIndx = out
  //   send_buffer[ offset + 0 ] = (Real) partIDs[pIndx];
  //   send_buffer[ offset + 1 ] = mass[ pIndx ];
  //   send_buffer[ offset + 2 ] = pos_x[ pIndx ];
  //   send_buffer[ offset + 3 ] = pos_y[ pIndx ];
  //   send_buffer[ offset + 4 ] = pos_z[ pIndx ];
  //   send_buffer[ offset + 5 ] = vel_x[ pIndx ];
  //   send_buffer[ offset + 6 ] = vel_y[ pIndx ];
  //   send_buffer[ offset + 7 ] = vel_z[ pIndx ];
  //   send_buffer[buffer_start] += 1;
  //   offset += N_DATA_PER_PARTICLE_TRANSFER;
  // }
  // send_buffer[buffer_start] = 10;

}


void Particles_3D::Add_Particle_To_Buffer( Real *buffer, int buffer_start, Real pId, Real pMass,
                            Real pPos_x, Real pPos_y, Real pPos_z,
                            Real pVel_x, Real pVel_y, Real pVel_z ){
  int n_in_buffer, offset_out;
  n_in_buffer = real_to_int( buffer[buffer_start] );
  offset_out = buffer_start + N_HEADER_PARTICLES_TRANSFER + n_in_buffer*N_DATA_PER_PARTICLE_TRANSFER;
  buffer[offset_out + 0] = pId;
  buffer[offset_out + 1] = pMass;
  buffer[offset_out + 2] = pPos_x;
  buffer[offset_out + 3] = pPos_y;
  buffer[offset_out + 4] = pPos_z;
  buffer[offset_out + 5] = pVel_x;
  buffer[offset_out + 6] = pVel_y;
  buffer[offset_out + 7] = pVel_z;
  buffer[buffer_start] += 1;
}

void Particles_3D::Add_Particle_To_Vectors( Real pId, Real pMass,
                            Real pPos_x, Real pPos_y, Real pPos_z,
                            Real pVel_x, Real pVel_y, Real pVel_z ){
  partIDs.push_back( real_to_part_int(pId) );
  mass.push_back( pMass );
  pos_x.push_back( pPos_x );
  pos_y.push_back( pPos_y );
  pos_z.push_back( pPos_z );
  vel_x.push_back( pVel_x );
  vel_y.push_back( pVel_y );
  vel_z.push_back( pVel_z );
  grav_x.push_back(0);
  grav_y.push_back(0);
  grav_z.push_back(0);

  n_local += 1;
}

void Particles_3D::Unload_Particles_from_Buffer( int direction, int side, int buffer_start, Real *recv_buffer, Real *send_buffer_y0, Real *send_buffer_y1, Real *send_buffer_z0, Real *send_buffer_z1, Real buffer_start_y, Real buffer_start_z  ){

  int n_recv;
  int offset_buff, indx, offset_out, n_in_out;
  part_int_t pId;
  Real pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z;

  n_recv = ( int ) recv_buffer[buffer_start];
  int n_added = 0;
  offset_buff = buffer_start + N_HEADER_PARTICLES_TRANSFER;
  for ( indx = 0; indx<n_recv; indx++ ){
    pId    = recv_buffer[ offset_buff + 0 ];
    pMass  = recv_buffer[ offset_buff + 1 ];
    pPos_x = recv_buffer[ offset_buff + 2 ];
    pPos_y = recv_buffer[ offset_buff + 3 ];
    pPos_z = recv_buffer[ offset_buff + 4 ];
    pVel_x = recv_buffer[ offset_buff + 5 ];
    pVel_y = recv_buffer[ offset_buff + 6 ];
    pVel_z = recv_buffer[ offset_buff + 7 ];
    offset_buff += N_DATA_PER_PARTICLE_TRANSFER;
    if ( pPos_x <  G.domainMin_x ) pPos_x += ( G.domainMax_x - G.domainMin_x );
    if ( pPos_x >= G.domainMax_x ) pPos_x -= ( G.domainMax_x - G.domainMin_x );
    if ( ( pPos_x < G.xMin ) || ( pPos_x >= G.xMax )  ){
      std::cout << "ERROR Particle Transfer out of X domain" << std::endl;
      std::cout << " posX: " << pPos_x << "velX: " << pVel_x << std::endl;
      std::cout << " posY: " << pPos_y << "velY: " << pVel_y << std::endl;
      std::cout << " posZ: " << pPos_z << "velZ: " << pVel_z << std::endl;
      continue;
    }

    if (  direction == 1 ){
      if ( pPos_y <  G.domainMin_y ) pPos_y += ( G.domainMax_y - G.domainMin_y );
      if ( pPos_y >= G.domainMax_y ) pPos_y -= ( G.domainMax_y - G.domainMin_y );
    }
    if ( (direction==1 || direction==2) && (( pPos_y < G.yMin ) || ( pPos_y >= G.yMax ))  ){
      std::cout << "ERROR Particle Transfer out of Y domain" << std::endl;
      std::cout << " posX: " << pPos_x << "velX: " << pVel_x << std::endl;
      std::cout << " posY: " << pPos_y << "velY: " << pVel_y << std::endl;
      std::cout << " posZ: " << pPos_z << "velZ: " << pVel_z << std::endl;
      continue;
    }

    if (direction  != 1 ){
      if ( pPos_y < G.yMin ){
        // std::cout << "  Adding to Y0 " << std::endl;
        Add_Particle_To_Buffer( send_buffer_y0, buffer_start_y, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
        continue;
      }
      if ( pPos_y >= G.yMax ){
        // std::cout << "  Adding to Y1 " << std::endl;
        Add_Particle_To_Buffer( send_buffer_y1, buffer_start_y, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
        continue;
      }
    }

    if (  direction == 2 ){
      if ( pPos_z <  G.domainMin_z ) pPos_z += ( G.domainMax_z - G.domainMin_z );
      if ( pPos_z >= G.domainMax_z ) pPos_z -= ( G.domainMax_z - G.domainMin_z );
    }
    if ( (direction==2) && (( pPos_z < G.zMin ) || ( pPos_z >= G.zMax ))  ){
      std::cout << "ERROR Particle Transfer out of Z domain" << std::endl;
      std::cout << " posX: " << pPos_x << "velX: " << pVel_x << std::endl;
      std::cout << " posY: " << pPos_y << "velY: " << pVel_y << std::endl;
      std::cout << " posZ: " << pPos_z << "velZ: " << pVel_z << std::endl;
      continue;
    }
  //
    if (direction  !=2 ){
      if ( pPos_z < G.zMin ){
        // std::cout << "  Adding to Z0 " << std::endl;
        Add_Particle_To_Buffer( send_buffer_z0, buffer_start_z, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
        continue;
      }
      if ( pPos_z >= G.zMax ){
        // std::cout << "  Adding to Z1 " << std::endl;
        Add_Particle_To_Buffer( send_buffer_z1, buffer_start_z, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
        continue;
      }
    }
    n_added += 1;
    Add_Particle_To_Vectors( pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
  }
  // std::cout << "N_recv: " << n_recv << "  N_added: " << n_added << std::endl;
}

void Particles_3D::Remove_Transfered_Particles( void ){
  // part_int_t n_delete = 0;
  // n_delete += out_indxs_vec_x0.size();
  // n_delete += out_indxs_vec_x1.size();
  // n_delete += out_indxs_vec_y0.size();
  // n_delete += out_indxs_vec_y1.size();
  // n_delete += out_indxs_vec_z0.size();
  // n_delete += out_indxs_vec_z1.size();
  //
  // int_vector_t delete_indxs_vec;
  // delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_x0.begin(), out_indxs_vec_x0.end() );
  // delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_x1.begin(), out_indxs_vec_x1.end() );
  // delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_y0.begin(), out_indxs_vec_y0.end() );
  // delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_y1.begin(), out_indxs_vec_y1.end() );
  // delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_z0.begin(), out_indxs_vec_z0.end() );
  // delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_z1.begin(), out_indxs_vec_z1.end() );
  //
  // out_indxs_vec_x0.clear();
  // out_indxs_vec_x1.clear();
  // out_indxs_vec_y0.clear();
  // out_indxs_vec_y1.clear();
  // out_indxs_vec_z0.clear();
  // out_indxs_vec_z1.clear();
  //
  //
  // std::sort(delete_indxs_vec.begin(), delete_indxs_vec.end());
  //
  // part_int_t indx, pIndx;
  // for ( indx=0; indx<n_delete; indx++ ){
  //   pIndx = delete_indxs_vec.back();
  //   Remove_ID( pIndx, partIDs );
  //   Remove_Real( pIndx, mass );
  //   Remove_Real( pIndx, pos_x );
  //   Remove_Real( pIndx, pos_y );
  //   Remove_Real( pIndx, pos_z );
  //   Remove_Real( pIndx, vel_x );
  //   Remove_Real( pIndx, vel_y );
  //   Remove_Real( pIndx, vel_z );
  //   Remove_Real( pIndx, grav_x );
  //   Remove_Real( pIndx, grav_y );
  //   Remove_Real( pIndx, grav_z );
  //   delete_indxs_vec.pop_back();
  //   n_local -= 1;
  // }
  // if ( delete_indxs_vec.size() != 0 ) std::cout << "ERROR: Deleting Transfered Particles " << std::endl;

  int n_in_out_vectors, n_in_vectors;
  n_in_vectors = ( partIDs.size() + mass.size() + pos_x.size() + pos_y.size() + pos_z.size() + vel_x.size() + vel_y.size() + vel_z.size() ) / N_DATA_PER_PARTICLE_TRANSFER;
  if ( n_in_vectors != n_local ) std::cout << "ERROR PARTICLES TRANSFER: DATA IN VECTORS DIFFERENT FROM N_LOCAL###########" << std::endl;
  n_in_out_vectors = out_indxs_vec_x0.size() + out_indxs_vec_x1.size() + out_indxs_vec_y0.size() + out_indxs_vec_y1.size() + out_indxs_vec_z0.size() + out_indxs_vec_z1.size();
  if ( n_in_out_vectors != 0 ) std::cout << "ERROR PARTICLES TRANSFER: OUPTUT VECTORS NOT EMPTY###########" << std::endl;
  // for ( int i=0; i<nproc; i++ ){
    // MPI_Barrier(world);
  // }
  // MPI_Barrier(world);
  // std::cout << "   Removed: "<< n_delete << std::endl;

}


#endif
#endif
