EXEC   = cholla

OPTIMIZE =  -O3

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

ifdef SELF_GRAVITY
PFFTINC = -I/home/bruno/apps/pfft-1.0.8-alpha/include
FFTWINC = -I/home/bruno/apps/fftw-3.3.5/include
PFFTLIB = -L/home/bruno/apps/pfft-1.0.8-alpha/lib
FFTWLIB = -L/home/bruno/apps/fftw-3.3.5/lib
PFFT_INC = $(PFFTINC) $(FFTWINC)
PFFT_LIBS = $(PFFTLIB) $(FFTWLIB)
PFFTFLAGS = $(PFFT_INC) -lpfft -lfftw3_mpi -lfftw3_omp -lfftw3 -fopenmp
endif

# Turn on particles
# PARTICLES = -DPARTICLES
# LONG_IDS = -DLONG_IDS
#
# # Turn on cosmology
# COSMOLOGY = -DCOSMOLOGY
# COSMO_UNITS = -DCOSMO_UNITS
#
# OMP_OVERLAP_HYDRO = -DOMP_OVERLAP_HYDRO

#PFFT_OMP = -DPFFT_OMP
ifdef PFFT_OMP
PFFTFLAGS += -lpfft_omp
endif




ifdef CUDA
CUDA_INCLUDE = -I/usr/local/cuda-9.0/include/
CUDA_LIBS = -L/usr/local/cuda-9.0/lib64/ -lcuda -lcudart
endif

ifeq ($(OUTPUT),-DHDF5)
HDF5_INCLUDE = -I/usr/include/hdf5/serial/
HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5
endif

INCL   = -I./ $(HDF5_INCLUDE)
NVINCL = $(INCL) $(CUDA_INCLUDE)
NVLIBS = $(CUDA_LIBS)
LIBS   = -lm $(HDF5_LIBS)



FLAGS = $(CUDA) $(PRECISION) $(OUTPUT) $(RECONSTRUCTION) $(SOLVER) $(INTEGRATOR) $(H_CORRECTION) $(COOLING) $(DUAL_ENERGY) $(SELF_GRAVITY) $(PARTICLES) $(LONG_IDS) $(COSMOLOGY) $(COSMO_UNITS) $(OMP_OVERLAP_HYDRO)
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
