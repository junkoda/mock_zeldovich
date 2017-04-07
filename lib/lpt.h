#ifndef LPT_H
#define LPT_H 1

#include "mem.h"
#include "particle.h"
#include "power.h"
#include "fft.h"

void lpt_init(const int nc, const double boxsize, Mem* mem);
void lpt_free();

void lpt_set_displacements(const unsigned long seed, PowerSpectrum* const ps,
			   const double a, Particles* particles);
FFT* lpt_generate_phi(const unsigned long seed, PowerSpectrum* const);
void lpt_set_offset(Float offset_);

#endif
