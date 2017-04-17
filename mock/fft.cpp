//
// Interface for FFTW
//

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <fftw3-mpi.h>
#include "config.h"
#include "comm.h"
#include "mem.h"
#include "msg.h"
#include "util.h"
#include "fft.h"

using namespace std;


FFT::FFT(const char name[], const int nc_, Mem* mem, const bool transposed) :
  nc(nc_), mode(fft_mode_unknown), own_mem(nullptr)
{
  // Allocates memory for FFT real and Fourier space and initilise fftw_plans
  assert(nc > 0);

  msg_printf(msg_verbose, "Setting up FFT %s with FFTW_MEASURE\n", name);
  
  if(transposed) {
    ncomplex= FFTW(mpi_local_size_3d_transposed)(nc, nc, nc/2+1, MPI_COMM_WORLD,
	                 &local_nx, &local_ix0,
			 &local_nky, &local_iky0);
  }
  else {
    ncomplex= FFTW(mpi_local_size_3d)(nc, nc, nc/2+1, MPI_COMM_WORLD,
			    &local_nx, &local_ix0);
    local_nky= local_iky0= 0;
  }

  size_t size= sizeof(complex_t)*ncomplex;
  assert(local_nx >= 0); assert(local_ix0 >= 0);
  
    
  if(mem == 0) {
    mem= new Mem(name, size);
    own_mem= mem;
  }
  
  void* buf= mem->use_remaining(size);
  // Call mem_use_from_zero(mem, 0) before this to use mem from the beginning.

  fx= (Float*) buf; fk= (complex_t*) buf;

  unsigned flag= 0;
  if(transposed) flag= FFTW_MPI_TRANSPOSED_OUT;
  forward_plan= FFTW(mpi_plan_dft_r2c_3d)(nc, nc, nc, fx, fk,
					  MPI_COMM_WORLD, FFTW_MEASURE | flag);

  unsigned flag_inv= 0;
  if(transposed) {
    flag_inv= FFTW_MPI_TRANSPOSED_IN;
    msg_printf(msg_debug, "FFTW transposed in/out\n");
  }
  
  inverse_plan= FFTW(mpi_plan_dft_c2r_3d)(nc, nc, nc, fk, fx,
                                     MPI_COMM_WORLD, FFTW_MEASURE | flag_inv);
}


FFT::~FFT()
{
  if(comm_status() == comm_parallel) {
    FFTW(destroy_plan)(forward_plan);
    FFTW(destroy_plan)(inverse_plan);
  }

  if(own_mem)
    delete own_mem;
}


void FFT::execute_forward()
{
  if(mode != fft_mode_x) {
    msg_printf(msg_warn,
	       "Warning: FFT %s mode is %d not %d for execute_forward\n",
	       name, mode, fft_mode_x);
    throw FFTError();
  }
  FFTW(mpi_execute_dft_r2c)(forward_plan, fx, fk);
}

void FFT::execute_inverse()
{
  if(mode != fft_mode_k) {
    msg_printf(msg_warn,
	       "Warning: FFT %s mode is %d not %d for execute_inverse\n",
	       name, mode, fft_mode_x);
    throw FFTError();
  }
  FFTW(mpi_execute_dft_c2r)(inverse_plan, fk, fx);
}


size_t fft_mem_size(const int nc, const int transposed)
{
  // return the memory size necessary for the 3D FFT
  ptrdiff_t local_nx, local_ix0, local_nky, local_iky0;

  ptrdiff_t n= 0;
  if(transposed)
    n= FFTW(mpi_local_size_3d_transposed)(nc, nc, nc/2+1, MPI_COMM_WORLD,
	           &local_nx, &local_ix0, &local_nky, &local_iky0);
  else
    n= FFTW(mpi_local_size_3d)(nc, nc, nc/2+1, MPI_COMM_WORLD,
				      &local_nx, &local_ix0);

  return size_align(sizeof(complex_t)*n);
}


size_t fft_local_nx(const int nc)
{
  ptrdiff_t local_nx, local_ix0; 
  FFTW(mpi_local_size_3d)(nc, nc, nc/2+1, MPI_COMM_WORLD, &local_nx, &local_ix0);

  return local_nx;
}


void fft_finalise()
{
  if(comm_status() == comm_parallel)
    FFTW(mpi_cleanup)();
}

void* fft_malloc(size_t size)
{
  return FFTW(malloc)(size);
}

// Quotes
// "it is probably better for you to simply create multiple plans
//  (creating a new plan is quick once one exists for a given size)
// -- FFTW3 Manual for version 3.3.3 section 4.6

// "To prevent memory leaks, you must still call fftw_destroy_plan
//  before executing fftw_cleanup."
