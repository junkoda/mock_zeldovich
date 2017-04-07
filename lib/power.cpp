//
// Power spectrum for initial condition
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <gsl/gsl_spline.h>
#include "comm.h"
#include "msg.h"
#include "power.h"
#include "error.h"

PowerSpectrum::PowerSpectrum(const char filename[]) :
  log_k(0), interp_(0), acc_(0)
{
  if(comm_this_node() == 0)
    read_file_(filename);

  comm_bcast_int(&n, 1);

  if(comm_this_node() != 0) {
    log_k= (double*) malloc(sizeof(double)*2*n);
    log_P= log_k + n;
  }

  comm_bcast_double(log_k, 2*n);

  
  interp_= gsl_interp_alloc(gsl_interp_cspline, n);
  acc_= gsl_interp_accel_alloc();

  const int n_required= (int) gsl_interp_min_size(interp_);
  if(n < n_required)
    msg_printf(msg_fatal, "Error: Not enough power spectrum data points for cubic spline; %d data points < %d required\n", n, n_required);

  gsl_interp_init(interp_, log_k, log_P, n);
}


PowerSpectrum::~PowerSpectrum()
{
  if(interp_) gsl_interp_free(interp_);
  if(acc_) gsl_interp_accel_free(acc_);
  if(log_k) free(log_k);
}


double PowerSpectrum::P(const double k) const
{
  double logP= gsl_interp_eval(interp_, log_k, log_P, log(k), acc_);
  return exp(logP);
}


void PowerSpectrum::read_file_(const char filename[])
{
  // Input:  filename of power spectrum (space separated ascii file)
  // Output: Arrays of log(k) and log(P) as this->log_k, this->log_P
  
  // File format, requirements, and assumptions
  //   - #: comment if first character is #
  //   - k P -- space spearated. OK to have 3 or more columns (neglected).
  //   - k must be in increasing order [1/h Mpc]
  //   - one line must be less than 128 characters, including \n.
  //   - no requirements for the number of lines.
  
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal, "Error: Unable to open input power spectrum file: %s\n",filename);
    throw IOError();
  }

  int nalloc= 1000;
    // Initial guess of number of data lines. Best if it is the number of the
    // power spectrum data entries, but works for any number.
  double* buf= (double*) malloc(sizeof(double)*2*nalloc);

  char line[128];
  int nlines= 0;
  double k, P;
  double k_prev= 0.0;

  // Read lines and push to buf as k1,P1,k2,P2, ...
  // Doubles the length of buf when the length of the array is not enough
  while(fgets(line, 127, fp)) {
    if(nlines == nalloc) {
      msg_printf(msg_debug, "reallocating power spectrum table %d -> %d\n",
		 nalloc, 2*nalloc);
      nalloc *= 2;
      buf= (double*) realloc(buf, sizeof(double)*2*nalloc); assert(buf);
    }
    
    if(line[0] == '#')
      continue;

    int ret= sscanf(line, "%lg %lg", &k, &P);
    if(ret != 2) {
      msg_printf(msg_warn,
		 "Warning: Unable to understand a line in the power spectrum "
		 "file; following data are ignored: %s", line);
      break;
    }

    if(k <= 0.0) {
      msg_printf(msg_warn, "Skipped power k= %f\n", k);
      continue;
    }
      
    if(k < k_prev) {
      msg_printf(msg_fatal,
		 "Error: wavenumber k in the power spectrum file must be sorted"
		 " in increasing order. %dth data k=%e > previous k= %e\n",
		 nlines, k_prev, k);
      throw IOError();
    }

    buf[2*nlines    ]= k;
    buf[2*nlines + 1]= P;

    k_prev= k;
    nlines++;
  }

  int ret= fclose(fp); assert(ret == 0);
  
  msg_printf(msg_verbose, "Read %d pairs of k P(k) from %s\n", nlines,filename);

  // Allocate ps->log_k, ps->log_P and fill the arrays
  double* const v_logk= (double*) malloc(2*nlines*sizeof(double));
  assert(v_logk);
  double* const v_logP= v_logk + nlines;

  for(int j=0; j<nlines; j++) {
    v_logk[j]= log(buf[2*j    ]);
    v_logP[j]= log(buf[2*j + 1]);
  }
  free(buf);
  
  log_k = v_logk;
  log_P = v_logP;
  n     = nlines;
}

double PowerSpectrum::compute_sigma(const double R) const
{
  // Computes sigma (rms amplituede) smoothed on scale R
  // R: smoothing length [/h Mpc] (8 for sigma_8)
  // 1/(2*pi^2) \int P(k) W(k*R) k^2 dk
  
  const double fac= 1.0/(2.0*M_PI*M_PI);

  double k0= exp(log_k[0]);
  double f0= 0.0;
  
  double sigma2= 0.0;
  for(int i=0; i<n; i++) {
    double k= exp(log_k[i]);
    double x= k*R;
    double w= 3.0*(sin(x)-x*cos(x))/(x*x*x);
    double f1= exp(log_P[i])*k*k*w*w;
    
    sigma2 += 0.5*(f0 + f1)*(k - k0);

    k0= k;
    f0= f1;
  }

  return sqrt(fac*sigma2);
}
