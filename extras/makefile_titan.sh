module switch PrgEnv-pgi/5.2.82 PrgEnv-gnu
module load cray-hdf5
module load cudatoolkit
module load fftw

EXEC   = cholla

OPTIMIZE =  -O2

DIR = ./src
CFILES = $(wildcard $(DIR)/*.c)
CPPFILES = $(wildcard $(DIR)/*.cpp)
CUDAFILES = $(wildcard $(DIR)/*.cu)


DIR_GRAV = ./src/gravity
CFILES_GRAV = $(wildcard $(DIR_GRAV)/*.c)
CPPFILES_GRAV = $(wildcard $(DIR_GRAV)/*.cpp)
CUDAFILES_GRAV = $(wildcard $(DIR_GRAV)/*.cu)

DIR_PART = ./src/particles
CFILES_PART = $(wildcard $(DIR_PART)/*.c)
CPPFILES_PART = $(wildcard $(DIR_PART)/*.cpp)
CUDAFILES_PART = $(wildcard $(DIR_PART)/*.cu)

DIR_COSMO = ./src/cosmology
CFILES_COSMO = $(wildcard $(DIR_COSMO)/*.c)
CPPFILES_COSMO = $(wildcard $(DIR_COSMO)/*.cpp)
CUDAFILES_COSMO = $(wildcard $(DIR_COSMO)/*.cu)


OBJS   = $(subst .c,.o,$(CFILES)) $(subst .cpp,.o,$(CPPFILES)) $(subst .cu,.o,$(CUDAFILES))  $(subst .c,.o,$(CFILES_GRAV)) $(subst .cpp,.o,$(CPPFILES_GRAV)) $(subst .cu,.o,$(CUDAFILES_GRAV)) $(subst .c,.o,$(CFILES_PART)) $(subst .cpp,.o,$(CPPFILES_PART)) $(subst .cu,.o,$(CUDAFILES_PART)) $(subst .c,.o,$(CFILES_COSMO)) $(subst .cpp,.o,$(CPPFILES_COSMO)) $(subst .cu,.o,$(CUDAFILES_COSMO))
COBJS   = $(subst .c,.o,$(CFILES)) $(subst .c,.o,$(CFILES_GRAV)) $(subst .c,.o,$(CFILES_PART)) $(subst .c,.o,$(CFILES_COSMO))
CPPOBJS   = $(subst .cpp,.o,$(CPPFILES)) $(subst .cpp,.o,$(CPPFILES_GRAV)) $(subst .cpp,.o,$(CPPFILES_PART)) $(subst .cpp,.o,$(CFILES_COSMO))
CUOBJS   = $(subst .cu,.o,$(CUDAFILES)) $(subst .cu,.o,$(CUDAFILES_GRAV)) $(subst .cu,.o,$(CUDAFILES_PART)) $(subst .cu,.o,$(CFILES_COSMO))

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
#RECONSTRUCTION = -DPPMP
RECONSTRUCTION = -DPPMC

#SOLVER = -DEXACT
#SOLVER = -DROE
SOLVER = -DHLLC

#INTEGRATOR = -DCTU
INTEGRATOR = -DVL

H_CORRECTION = -DH_CORRECTION
#COOLING = -DCOOLING_GPU

DUAL_ENERGY = -DDE


#Turn on SELF_GRAVITY
GRAVITY = -DGRAVITY
#GRAVITY_RIEMANN = -DGRAVITY_RIEMANN
#POTENTIAL = -DPOTENTIAL_CUFFT
#POTENTIAL = -DPOTENTIAL_FFTW
POTENTIAL = -DPOTENTIAL_PFFT
#OUTPUT_POTENTIAL = -DOUTPUT_POTENTIAL
#OUTPUT_GRAVITY_DENSITY = -DOUTPUT_GRAVITY_DENSITY

#Turn on Particle-Mesh scheme
PARTICLES = -DPARTICLES
PARTICLE_INT = -DLONG_INTS
# PARTICLES_KDK = -DPARTICLES_KDK
# PARTICLES_DKD = -DPARTICLES_DKD
PARTICLES_OMP = -DPARTICLES_OMP
N_OMP_PARTICLE_THREADS = -DN_OMP_PARTICLE_THREADS=16
# PARTICLES_CUDA = -DPARTICLES_CUDA
# ONLY_PM = -DONLY_PM
# PECULIAR_VEL = -DPECULIAR_VEL
SINGLE_PARTICLE_MASS = -DSINGLE_PARTICLE_MASS

#Turn On Cosmological simulation0
COSMOLOGY = -DCOSMOLOGY


ifdef CUDA
CUDA_INCL = -I$(CUDATOOLKIT_HOME)/include/
CUDA_LIBS = -L$(CUDATOOLKIT_HOME)/lib64/ -lcuda -lcudart
endif

ifeq ($(OUTPUT),-DHDF5)
HDF5_INCL = -I$(HDF5_DIR)/include
HDF5_LIBS = -L$(HDF5_DIR)/lib -lhdf5
endif

INCL   = -I./ $(HDF5_INCL)
NVINCL = $(INCL) $(CUDA_INCL)
LIBS   = -lm $(HDF5_LIBS) $(CUDA_LIBS)

FFTW_INCL = -I/opt/cray/fftw/3.3.4.11/interlagos/include
FFTW_LIBS = -L/opt/cray/fftw/3.3.4.11/interlagos/lib -lfftw3

ifeq ($(POTENTIAL),-DPOTENTIAL_PFFT)
PFFT_INCL = -I/ccs/home/bvilasen/apps/PFFT/pfft-git/include
PFFT_LIBS = -L/ccs/home/bvilasen/apps/PFFT/pfft-git/lib  -lpfft  -lfftw3_mpi -lfftw3
INCL += $(FFTW_INCL) $(PFFT_INCL)
LIBS += $(FFTW_LIBS) $(PFFT_LIBS)
endif

ifdef PARTICLES_OMP
OMP_FLAGS = -fopenmp
LIBS += -fopenmp
endif

FLAGS = $(CUDA) $(PRECISION) $(OUTPUT) $(RECONSTRUCTION) $(SOLVER) $(INTEGRATOR) $(H_CORRECTION) $(COOLING) $(DUAL_ENERGY) $(GRAVITY) $(GRAVITY_RIEMANN) $(POTENTIAL) $(OUTPUT_POTENTIAL) $(OUTPUT_GRAVITY_DENSITY) $(PARTICLES)  $(SINGLE_PARTICLE_MASS) $(PARTICLE_INT) $(PARTICLES_KDK) $(PECULIAR_VEL) $(PARTICLES_OMP) $(PARTICLES_CUDA) $(N_OMP_PARTICLE_THREADS) $(ONLY_PM) $(COSMOLOGY)#-DSTATIC_GRAV #-DDE -DSCALAR -DSLICES -DPROJECTION -DROTATED_PROJECTION
CFLAGS 	  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS) $(OMP_FLAGS)
CXXFLAGS  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS) $(OMP_FLAGS)
NVCCFLAGS = $(FLAGS) -m64 -arch=compute_35 -fmad=false -ccbin=$(CC) --expt-relaxed-constexpr


%.o:	%.c
		$(CC) $(CFLAGS)  $(INCL)  -c $< -o $@

%.o:	%.cpp
		$(CXX) $(CXXFLAGS)  $(INCL) $(CUDA_INCL) -c $< -o $@

%.o:	%.cu
		$(NVCC) $(NVCCFLAGS) --device-c $(NVINCL)  -c $< -o $@

$(EXEC): $(OBJS) src/gpuCode.o
	 	 $(CXX) $(OBJS) src/gpuCode.o $(LIBS) -o $(EXEC)

src/gpuCode.o:	$(CUOBJS)
		$(NVCC)  -dlink $(CUOBJS) -o src/gpuCode.o



.PHONY : clean

clean:
	 rm -f $(OBJS) src/gpuCode.o $(EXEC)
