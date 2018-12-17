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

DIR_COOL = ./src/cooling
CFILES_COOL = $(wildcard $(DIR_COOL)/*.c)
CPPFILES_COOL = $(wildcard $(DIR_COOL)/*.cpp)
CUDAFILES_COOL = $(wildcard $(DIR_COOL)/*.cu)


OBJS   = $(subst .c,.o,$(CFILES)) $(subst .cpp,.o,$(CPPFILES)) $(subst .cu,.o,$(CUDAFILES))  $(subst .c,.o,$(CFILES_GRAV)) $(subst .cpp,.o,$(CPPFILES_GRAV)) $(subst .cu,.o,$(CUDAFILES_GRAV)) $(subst .c,.o,$(CFILES_PART)) $(subst .cpp,.o,$(CPPFILES_PART)) $(subst .cu,.o,$(CUDAFILES_PART)) $(subst .c,.o,$(CFILES_COSMO)) $(subst .cpp,.o,$(CPPFILES_COSMO)) $(subst .cu,.o,$(CUDAFILES_COSMO)) $(subst .c,.o,$(CFILES_COOL)) $(subst .cpp,.o,$(CPPFILES_COOL)) $(subst .cu,.o,$(CUDAFILES_COOL))
COBJS   = $(subst .c,.o,$(CFILES)) $(subst .c,.o,$(CFILES_GRAV)) $(subst .c,.o,$(CFILES_PART)) $(subst .c,.o,$(CFILES_COSMO)) $(subst .c,.o,$(CFILES_COOL))
CPPOBJS   = $(subst .cpp,.o,$(CPPFILES)) $(subst .cpp,.o,$(CPPFILES_GRAV)) $(subst .cpp,.o,$(CPPFILES_PART)) $(subst .cpp,.o,$(CPPFILES_COSMO)) $(subst .cpp,.o,$(CPPFILES_COOL))
CUOBJS   = $(subst .cu,.o,$(CUDAFILES)) $(subst .cu,.o,$(CUDAFILES_GRAV)) $(subst .cu,.o,$(CUDAFILES_PART)) $(subst .cu,.o,$(CUDAFILES_COSMO)) $(subst .cu,.o,$(CUDAFILES_COOL))

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
#RECONSTRUCTION = -DPPMP
RECONSTRUCTION = -DPPMC

#SOLVER = -DEXACT
#SOLVER = -DROE
SOLVER = -DHLLC

#INTEGRATOR = -DCTU
INTEGRATOR = -DVL

H_CORRECTION = -DH_CORRECTION
#COOLING = #-DCOOLING_GPU -DCLOUDY_COOL

DUAL_ENERGY = -DDE


#Turn on SELF_GRAVITY
GRAVITY = -DGRAVITY
#GRAVITY_RIEMANN = -DGRAVITY_RIEMANN
#POTENTIAL = -DPOTENTIAL_CUFFT
#POTENTIAL = -DPOTENTIAL_FFTW
POTENTIAL = -DPOTENTIAL_PFFT
#OUTPUT_POTENTIAL = -DOUTPUT_POTENTIAL
#OUTPUT_GRAVITY_DENSITY = -DOUTPUT_GRAVITY_DENSITY
GRAVITY_OMP = -DGRAVITY_OMP
N_OMP_GRAVITY_THREADS = -DN_OMP_GRAVITY_THREADS=5
#GRAVITY_CORRECTOR = -DGRAVITY_CORRECTOR
#REVERT_STEP = -DREVERT_STEP

#Turn on Particle-Mesh scheme
PARTICLES = -DPARTICLES
PARTICLE_INT = -DLONG_INTS
#ONLY_PM = -DONLY_PM
SINGLE_PARTICLE_MASS = -DSINGLE_PARTICLE_MASS
PARTICLE_IDS = -DPARTICLE_IDS

#Turn On Cosmological simulation0
COSMOLOGY = -DCOSMOLOGY

#Turn on scalar fields
#SCALAR = -DSCALAR

#EXTRA_FIELD = -DEXTRA_FIELD

#COOLING = -DCOOLING_GRACKLE
#GRACKLE_PRECISION = -DCONFIG_BFLOAT_8

CPU_TIME = -DCPU_TIME

ifdef CUDA
CUDA_INCL = -I/usr/local/cuda-9.0/targets/x86_64-linux/include/
CUDA_LIBS = -L/usr/local/cuda-9.0/targets/x86_64-linux/lib/ -lcuda -lcudart
ifeq ($(POTENTIAL),-DPOTENTIAL_CUFFT)
CUDA_LIBS += -lcufft
endif
endif
ifeq ($(OUTPUT),-DHDF5)
HDF5_INCL = -I/usr/include/hdf5/serial/
HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5
endif

INCL   = -I./ $(HDF5_INCL)
NVINCL = $(INCL) $(CUDA_INCL)
LIBS   = -lm $(HDF5_LIBS) $(CUDA_LIBS)

FFTW_INCL = -I/home/bruno/apps/fftw-3.3.5/include
FFTW_LIBS = -L/home/bruno/apps/fftw-3.3.5/lib -lfftw3
ifeq ($(POTENTIAL),-DPOTENTIAL_FFTW)
INCL += $(FFTW_INCL)
LIBS += $(FFTW_LIBS)
endif

ifeq ($(POTENTIAL),-DPOTENTIAL_PFFT)
PFFT_INCL = -I/home/bruno/apps/pfft-1.0.8-alpha/include
PFFT_LIBS = -L/home/bruno/apps/pfft-1.0.8-alpha/lib  -lpfft  -lfftw3_mpi -lfftw3
INCL += $(FFTW_INCL) $(PFFT_INCL)
LIBS += $(FFTW_LIBS) $(PFFT_LIBS)
endif

ifeq ($(COOLING),-DCOOLING_GRACKLE)
GRACKLE_INCL = -I/home/bruno/local/include
GRACKLE_LIBS = -L/home/bruno/local/lib -lgrackle
INCL += $(GRACKLE_INCL)
LIBS += $(GRACKLE_LIBS)
endif


ifdef GRAVITY_OMP
OMP_FLAGS = -fopenmp
LIBS += -fopenmp
endif



FLAGS = $(CUDA) $(PRECISION) $(CPU_TIME) $(OUTPUT) $(RECONSTRUCTION) $(SOLVER) $(INTEGRATOR) $(H_CORRECTION) $(COOLING) $(DUAL_ENERGY) $(GRAVITY) $(GRAVITY_RIEMANN) $(POTENTIAL) $(OUTPUT_POTENTIAL) $(OUTPUT_GRAVITY_DENSITY) $(PARTICLES)  $(SINGLE_PARTICLE_MASS) $(PARTICLE_INT) $(PARTICLE_IDS) $(PARTICLES_KDK) $(PECULIAR_VEL) $(GRAVITY_OMP) $(N_OMP_GRAVITY_THREADS) $(ONLY_PM) $(COSMOLOGY) $(GRACKLE_PRECISION) $(GRAVITY_CORRECTOR) $(SCALAR) $(REVERT_STEP) $(EXTRA_FIELD)#-DSTATIC_GRAV #-DDE -DSCALAR -DSLICES -DPROJECTION -DROTATED_PROJECTION
CFLAGS 	  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS) $(OMP_FLAGS)
CXXFLAGS  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS) $(OMP_FLAGS)
NVCCFLAGS = $(FLAGS) -fmad=false -ccbin=$(CC)


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
