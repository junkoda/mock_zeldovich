
# Set compiler
CXX := mpic++ -std=c++11
OPT := -DDOUBLEPRECISION

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
# OPENMP  = -fopenmp
# -openmp for Intel icpc, -fopenmp for gcc
# MacOS build-in LLVM doesn't support OpenMPx

# Define paths of FFTW3 & GSL libraries if necessary.
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl
HDF5P_DIR ?= # parallel HDF5 library; e.g., brew install hdf5 --with-mpi

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR) $(HDF5P_DIR)

INCLS     += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

export CXX OPT OPENMP INCLS LIBS

all: lib mock 

.PHONY: lib mock clean

lib:
	cd lib && $(MAKE) lib

mock:
	cd mock && $(MAKE) mock

clean:
	cd lib && $(MAKE) clean
	cd mock && $(MAKE) clean
