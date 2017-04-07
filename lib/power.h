#ifndef POWER_H
#define POWER_H 1

#include <gsl/gsl_spline.h>

class PowerSpectrum {
 public:
  PowerSpectrum(const char filename[]);
  ~PowerSpectrum();
  double P(const double k) const;
  double compute_sigma(const double R) const;
    
  int n;
  double* log_k;
  double* log_P;
  
 private:
  gsl_interp *interp_;
  gsl_interp_accel *acc_;

  void read_file_(const char filename[]);    
};

#endif
