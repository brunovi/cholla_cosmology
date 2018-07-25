#ifdef COSMOLOGY
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>
#ifdef HDF5
#include<hdf5.h>
#endif
#include"io.h"
#include"grid3D.h"

#include"io_cosmo.h"

#include <iostream>
#include <fstream>
using namespace std;


int Load_scale_outputs( struct parameters P ) {
  char filename_1[100];
  // create the filename to read from
  strcpy(filename_1, P.times_output);

  string line;
  // ifstream myfile_out ("outputs_z0_n25.txt");
  ifstream myfile_out ( filename_1 );


  Real a_value;

  scale_outputs.clear();
  n_outputs = 0;

  if (myfile_out.is_open())
  {
    while ( getline (myfile_out,line) )
    {
      a_value = atof( line.c_str() );
      scale_outputs.push_back( a_value );
      n_outputs += 1;
      // chprintf("%f\n", a_value);
    }
    myfile_out.close();

    chprintf(" Loaded %d scale outputs \n", n_outputs);
  }

  else chprintf("Unable to open cosmology outputs file\n");

  return 0;
}





#endif
