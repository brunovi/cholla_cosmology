EXEC   = cholla

OPTIMIZE =  -O2

DIR = ./src
CFILES = $(wildcard $(DIR)/*.c)
CPPFILES = $(wildcard $(DIR)/*.cpp)
CUDAFILES = $(wildcard $(DIR)/*.cu)


OBJS   = $(subst .c,.o,$(CFILES)) $(subst .cpp,.o,$(CPPFILES)) $(subst .cu,.o,$(CUDAFILES))
COBJS   = $(subst .c,.o,$(CFILES))
CPPOBJS   = $(subst .cpp,.o,$(CPPFILES))
CUOBJS   = $(subst .cu,.o,$(CUDAFILES))

#To use GPUs, CUDA must be turned on here
#Optional error checking can also be enabled
CUDA = -DCUDA #-DCUDA_ERROR_CHECK

#To use MPI, MPI_FLAGS must be set to -DMPI_CHOLLA
#otherwise gcc/g++ will be used for serial compilation
MPI_FLAGS =  -DMPI_CHOLLA

ifdef MPI_FLAGS
  CC	= mpicc
  CXX   = mpicxx

  #MPI_FLAGS += -DSLAB
  MPI_FLAGS += -DBLOCK

else
  CC	= gcc
  CXX   = g++
endif

#define the NVIDIA CUDA compiler
NVCC	= nvcc

.SUFFIXES : .c .cpp .cu .o

#PRECISION = -DPRECISION=1
PRECISION = -DPRECISION=2

#OUTPUT = -DBINARY
OUTPUT = -DHDF5

#RECONSTRUCTION = -DPCM
#RECONSTRUCTION = -DPLMP
RECONSTRUCTION = -DPLMC
#RECONSTRUCTION = -DPPMP
#RECONSTRUCTION = -DPPMC

#SOLVER = -DEXACT
#SOLVER = -DROE
SOLVER = -DHLLC

#INTEGRATOR = -DCTU
INTEGRATOR = -DVL

COOLING = #-DCOOLING_GPU -DCLOUDY_COOL


ifdef CUDA
CUDA_INCLUDE = -I$(CUDATOOLKIT_HOME)/include/
CUDA_LIBS = -L$(CUDATOOLKIT_HOME)/lib64/ -lcuda -lcudart
endif

ifeq ($(OUTPUT),-DHDF5)
HDF5_INCL = -I$(HDF5_DIR)/include
HDF5_LIBS = -L$(HDF5_DIR)/lib -lhdf5
endif

INCL   = -I./ $(HDF5_INCL)
NVINCL = $(INCL) $(CUDA_INCL)
LIBS   = -lm $(HDF5_LIBS) $(CUDA_LIBS)


FLAGS = $(CUDA) $(PRECISION) $(OUTPUT) $(RECONSTRUCTION) $(SOLVER) $(INTEGRATOR) $(COOLING) #-DSTATIC_GRAV #-DDE -DSCALAR -DSLICES -DPROJECTION -DROTATED_PROJECTION
CFLAGS 	  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS)
CXXFLAGS  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS)
NVCCFLAGS = $(FLAGS) -m64 -fmad=false -ccbin=$(CC) -arch=compute_35


%.o:	%.c
		$(CC) $(CFLAGS)  $(INCL)  -c $< -o $@

%.o:	%.cpp
		$(CXX) $(CXXFLAGS)  $(INCL) -c $< -o $@

%.o:	%.cu
		$(NVCC) $(NVCCFLAGS) --device-c $(NVINCL)  -c $< -o $@

src/gpuCode.o:	$(CUOBJS)
$(NVCC) -arch=compute_35 -dlink $(CUOBJS) -o src/gpuCode.o

$(EXEC): $(OBJS) src/gpuCode.o
	 	 $(CXX) $(OBJS) src/gpuCode.o $(LIBS) -o $(EXEC)




.PHONY : clean

clean:
	 rm -f $(OBJS) src/gpuCode.o $(EXEC)
