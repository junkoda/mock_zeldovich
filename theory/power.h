#ifndef POWER_H
#define POWER_H 1

#include <gsl/gsl_spline.h>

typedef struct {
  int n;
  double* k;
  double* P;
  double* log_k;
  double* log_P;           // logP= log_k + n
  gsl_interp *interp;
  gsl_interp_accel *acc;
  double k_min, k_max;
} PowerSpectrum;

PowerSpectrum* power_alloc(const char filename[], const double sigma_8= 0.0);
PowerSpectrum* power_alloc_3(const char filename[]);
void power_free(PowerSpectrum* ps);

//double power(PowerSpectrum const * const ps, const double k);
double power_sigma8(PowerSpectrum const * const ps);
double power_sigmav(PowerSpectrum const * const ps);
void power_rescale(PowerSpectrum * const ps, const double growth);


#endif
