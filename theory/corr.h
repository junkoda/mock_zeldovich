#ifndef CORR_H
#define CORR_H 1

#include <gsl/gsl_spline.h>
#include "power.h"

typedef struct {
  int n;
  double const * r;
  double* xi;
  gsl_interp *interp;
  gsl_interp_accel *acc;  
} CorrelationFunction;

void corr_init(double const * const r, const size_t n,
	       PowerSpectrum const * const);
CorrelationFunction* corr_alloc(double const * const r, const size_t n);

void corr_integrate_psi_hybrid();
void corr_print();

void corr_fourier_test();
/*
double corr_dd(const int i);
double corr_du1(const int i);
double corr_uu0(const int i);
double corr_uu2(const int i);

double corr_dd(const double r);
double corr_du1(const double r);
double corr_uu0(const double r);
double corr_uu2(const double r);

CorrelationFunction* corr_alloc(double const * const r, const size_t n);

CorrelationFunction const * get_xi_dd();
CorrelationFunction const * get_xi_du1();
CorrelationFunction const * get_xi_uu0();
CorrelationFunction const * get_xi_uu2();
*/

#endif
