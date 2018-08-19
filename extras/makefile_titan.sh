EXEC   = cholla

OPTIMIZE =  -O2

DIR = ./src
CFILES = $(wildcard $(DIR)/*.c)
CPPFILES = $(wildcard $(DIR)/*.cpp)
CUDAFILES = $(wildcard $(DIR)/*.cu)


OBJS   = $(subst .c,.o,$(CFILES)) $(subst .cpp,.o,$(CPPFILES)) $(subst .cu,.o,$(CUDAFILES))

#To use GPUs, CUDA must be turned on here
#Optional error checking can also be enabled
CUDA = -DCUDA #-DCUDA_ERROR_CHECK

#To use MPI, MPI_FLAGS must be set to -DMPI_CHOLLA
#otherwise gcc/g++ will be used for serial compilation
MPI_FLAGS =  -DMPI_CHOLLA

ifdef MPI_FLAGS
  CC	= cc
  CXX   = CC

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
#RECONSTRUCTION = -DPLMC
RECONSTRUCTION = -DPPMP
#RECONSTRUCTION = -DPPMC

#SOLVER = -DEXACT
#SOLVER = -DROE
SOLVER = -DHLLC

#INTEGRATOR = -DCTU
INTEGRATOR = -DVL

H_CORRECTION = -DH_CORRECTION
#COOLING = -DCOOLING_GPU


DUAL_ENERGY = -DDE


# #Turn on SELF_GRAVITY
# SELF_GRAVITY = -DSELF_GRAVITY
#
# ifdef SELF_GRAVITY
# PFFTINC = /ccs/home/bvilasen/apps/PFFT/pfft-git/include
# PFFTLIB = /ccs/home/bvilasen/apps/PFFT/pfft-git/lib
# FFTWINC = /opt/cray/fftw/3.3.4.11/interlagos/include
# FFTWLIB = /opt/cray/fftw/3.3.4.11/interlagos/lib
# PFFT_INC = -I$(PFFTINC) -I$(FFTWINC)
# PFFT_LIBS = -L$(PFFTLIB) -L$(FFTWLIB)
# PFFTFLAGS = $(PFFT_INC) -openmp -lpfft -lfftw3_mpi -lfftw3
# endif
#
# PFFT_OMP = -DPFFT_OMP
# ifdef PFFT_OMP
# PFFTFLAGS += -lpfft_omp
# endif
#
# # Turn on particles
# PARTICLES = -DPARTICLES
# LONG_IDS = -DLONG_IDS
#
# # Turn on cosmology
# COSMOLOGY = -DCOSMOLOGY
# COSMO_UNITS = -DCOSMO_UNITS
#
# OMP_OVERLAP_HYDRO = -DOMP_OVERLAP_HYDRO


ifdef CUDA
CUDA_INCLUDE = -I$(CUDATOOLKIT_HOME)/include/
CUDA_LIBS = -L$(CUDATOOLKIT_HOME)/lib64/ -lcuda -lcudart
endif

ifeq ($(OUTPUT),-DHDF5)
HDF5_INCLUDE = -I$(HDF5_DIR)/include
HDF5_LIBS = -L$(HDF5_DIR)/lib -lhdf5
endif

INCL   = -I./ $(HDF5_INCLUDE)
NVINCL = $(INCL) $(CUDA_INCLUDE)
NVLIBS = $(CUDA_LIBS)
LIBS   = -lm $(HDF5_LIBS)



FLAGS = $(CUDA) $(PRECISION) $(OUTPUT) $(RECONSTRUCTION) $(SOLVER) $(INTEGRATOR) $(H_CORRECTION) $(COOLING) $(DUAL_ENERGY) $(SELF_GRAVITY) $(PARTICLES) $(LONG_IDS) $(COSMOLOGY) $(COSMO_UNITS) $(OMP_OVERLAP_HYDRO) $(PFFT_OMP)
CFLAGS 	  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS) $(FFT_FLAGS) -m64
CXXFLAGS  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS) $(FFT_FLAGS) $(PFFTFLAGS) -m64
NVCCFLAGS = $(FLAGS) -m64 -arch=compute_35 -fmad=false -ccbin=$(CC)
LDFLAGS	  = -m64


%.o:	%.c
		$(CC) $(CFLAGS)  $(INCL)  -c $< -o $@

%.o:	%.cpp
		$(CXX) $(CXXFLAGS)  $(INCL)  -c $< -o $@

%.o:	%.cu
		$(NVCC) $(NVCCFLAGS)  $(INCL)  -c $< -o $@


$(EXEC): $(OBJS)
	 	 $(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) $(NVLIBS) $(PFFTFLAGS) $(PFFT_LIBS) -o $(EXEC) $(INCL)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)
