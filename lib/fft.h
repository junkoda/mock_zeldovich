#ifndef FFT_H
#define FFT_H 1

#include <stdbool.h>
#include "config.h"
#include "mem.h"

// FFTW() adds fftw_ or fftwf_ prefix depending on DOUBLEPRECISION
#ifdef DOUBLEPRECISION
#define FFTW(f) fftw_ ## f
#else
#define FFTW(f) fftwf_ ## f
#endif


enum FFTMode {fft_mode_unknown, fft_mode_x, fft_mode_k};

class FFT {
 public:
  FFT(const char name[], const int nc, Mem* mem, const bool transposed);
  ~FFT();
  void execute_forward();
  void execute_inverse();
  
  char*       name;
  int         nc;
  Float*    fx;
  complex_t*  fk;
  ptrdiff_t   local_nx, local_ix0;
  ptrdiff_t   local_nky, local_iky0;
  FFTMode     mode;
 private:
  FFTW(plan)  forward_plan, inverse_plan;
  ptrdiff_t   ncomplex;
  Mem*        own_mem; // Allocated memory soley for this FFT
};

class FFTError{};

size_t fft_mem_size(const int nc, const int transposed);
size_t fft_local_nx(const int nc);
  
void fft_finalize();
void* fft_malloc(size_t size);

#endif
