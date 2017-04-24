#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_spline.h>
#include "power.h"

PowerSpectrum* power_alloc(const char filename[], const double sigma_8)
{
  PowerSpectrum* ps= (PowerSpectrum*) malloc(sizeof(PowerSpectrum)); assert(ps);

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    fprintf(stderr, "Error: Unable to open input power spectrum file: %s\n",filename);
    abort();
  }
      

  int nalloc= 1000;
  double* buf= (double*) malloc(sizeof(double)*2*nalloc);

  char line[128];
  int nlines= 0;
  double k, P;
  double k_prev= 0.0, P_prev= 0.0, sigma8_sq= 0.0;

  const double Rth= 8.0; // 8/h Mpc smoothing for sigma_8
  const double factor= 1.0/(2.0*M_PI*M_PI);

  // Read lines and push to buf as k1,P1,k2,P2, ...
  // Doubles the length of buf when the length of the array is not enough
  while(fgets(line, 127, fp)) {
    if(nlines == nalloc) {
      fprintf(stderr, "Reallocating power spectrum table %d -> %d\n",
		 nalloc, 2*nalloc);
      nalloc *= 2;
      buf= (double*) realloc(buf, sizeof(double)*2*nalloc); assert(buf);
    }
    
    if(line[0] == '#')
      continue;
    else if(sscanf(line, "%lg %lg", &k, &P) == 2) {
      if(k < k_prev) {
	fprintf(stderr, "Error: wavenumber k in the power spectrum file must be sorted in increasing order. %dth data k=%e > previous k= %e\n", nlines, k_prev, k);
	abort();
      }
      
      buf[2*nlines    ]= k;
      buf[2*nlines + 1]= P;

      // sigma_8 computation
      double x= k*Rth;
      double w= 3.0*(sin(x)-x*cos(x))/(x*x*x);

      sigma8_sq += 0.5*factor*(P*k*k + P_prev*k_prev*k_prev)*(w*w)*(k-k_prev);
      
      k_prev= k;
      P_prev= P;
      nlines++;
    }
    else {
      fprintf(stderr, "Warning: Unable to understand a line in the power spectrum file; following data are ignored: %s", line);
      break;
    }
  }

  int ret= fclose(fp); assert(ret == 0);
  
  fprintf(stderr, "Read %d pairs of k P(k) from %s\n", nlines, filename);

  // Check sigma8 (check skipped if sigma8_check = 0);
  const double sigma8_file= sqrt(sigma8_sq);

  double norm= 1.0;
  if(sigma_8 > 0.0) {
    fprintf(stderr, "normalize the power spectrum %e -> %e\n", sigma8_file, sigma_8);
    norm= (sigma_8/sigma8_file)*(sigma_8/sigma8_file);
  }
  else {
    fprintf(stderr, "sigma8= %f\n", sigma8_file);
  }
  
  // Allocate ps->log_k, ps->log_P and fill the arrays
  //double* const v_logk=
  ps->k   = (double*) malloc(4*nlines*sizeof(double)); assert(ps->k);
  ps->P   = ps->k + nlines;
  ps->log_k= ps->P + nlines;
  ps->log_P= ps->log_k + nlines;
    
  for(int j=0; j<nlines; j++) {
    ps->k[j] = buf[2*j];
    ps->P[j] = norm*buf[2*j + 1];
    ps->log_k[j]= log(buf[2*j]);
    ps->log_P[j]= log(norm*buf[2*j + 1]);
  }
  free(buf);
  
  ps->n= nlines;

  ps->interp= gsl_interp_alloc( gsl_interp_linear, ps->n);
 
  //ps->interp= gsl_interp_alloc(gsl_interp_cspline, ps->n);
  ps->acc= gsl_interp_accel_alloc();

  const int n_required= (int) gsl_interp_min_size(ps->interp);
  if(ps->n < n_required) {
    fprintf(stderr, "Error: Not enough power spectrum data points for cubic spline; %d data points < %d required\n", ps->n, n_required);
    abort();
  }

  gsl_interp_init(ps->interp, ps->log_k, ps->log_P, ps->n);

  ps->k_min= ps->k[0];
  ps->k_max= ps->k[nlines-1];

  return ps;
}

PowerSpectrum* power_alloc(const double k_min, const double k_max,
			   const double dk)
{
  const int n= round((k_max - k_min)/dk);
  
  PowerSpectrum* const ps=
    (PowerSpectrum*) malloc(sizeof(PowerSpectrum)); assert(ps);

  ps->k   = (double*) malloc(2*n*sizeof(double)); assert(ps->k);
  ps->P   = ps->k + n;
    
  for(int i=0; i<n; i++) {
    ps->k[i] = k_min + i*dk;
    ps->P[i] = 0.0;
  }
  
  ps->n= n;

  ps->k_min= ps->k[0];
  ps->k_max= ps->k[n-1];

  return ps;
}


void power_free(PowerSpectrum* const ps)
{
  gsl_interp_accel_free(ps->acc);
  gsl_interp_free(ps->interp);
  free(ps->k);
  free(ps);
}  


double power(PowerSpectrum const * const ps, const double k)
{
  if(k <= ps->k_min || k > ps->k_max)
    return 0.0;
  
  double log_P=
    gsl_interp_eval(ps->interp, ps->log_k, ps->log_P, log(k), ps->acc);

  return exp(log_P);
}


double power_sigma8(PowerSpectrum const * const ps)
{
  const int n= ps->n;
  const double Rth= 8.0; // 8/h Mpc smoothing for sigma_8
  const double factor= 1.0/(2.0*M_PI*M_PI);

  double f_prev= 0.0;
  double k_prev= ps->k[0];

  double sigma2= 0.0;
  
  for(int i=0; i<n; ++i) {
    double k= ps->k[i];
    double x= k*Rth;
    double w= 3.0*(sin(x)-x*cos(x))/(x*x*x);
    double f= ps->P[i]*k*k*w*w;
  
    sigma2 += 0.5*(f + f_prev)*(k-k_prev);

    k_prev= k;
    f_prev= f;
  }
  return sqrt(factor*sigma2);
}

double power_sigmav(PowerSpectrum const * const ps)
{
  const int n= ps->n;
  const double fac_uu= 1.0/(6.0*M_PI*M_PI);
  double sigma2= 0.0;
  
  for(int i=1; i<n; ++i) {
    double k1= ps->k[i-1];
    double k2= ps->k[i];
    
    sigma2 += 0.5*(ps->P[i-1] + ps->P[i])*(k2 - k1);
  }

  return sqrt(fac_uu*sigma2);
}

void power_rescale(PowerSpectrum * const ps, const double growth)
{
  const int n= ps->n;
  const double fac= growth*growth;
  
  for(int i=0; i<n; ++i) {
    ps->P[i] *= fac;
  }
}
