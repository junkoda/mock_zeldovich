//
// compute power spectrum from 3D delta(k)
//
#include <iostream>
#include <complex>
#include <vector>
#include <valarray>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <boost/program_options.hpp>

#include "msg.h"
#include "config.h"

using namespace std;
using namespace boost::program_options;


static vector<double> wcorr;

inline double w(const double tx)
{
  return tx == 0.0 ? 1.0 : sin(tx)/tx;
}

static void init_wcorr(const int nc, const int nmas)
{
  wcorr.resize(nc, 1.0);

  if(nmas == 0) return;

  for(int i=1; i<nc; ++i) {
    int k= i < nc/2 ? i : i - nc;
    
    wcorr[i]= pow(w(M_PI/nc*k), 2*nmas);
  }
}



void compute_power_spectrum(const char filename[],
		    complex_t* fk,
		    const int nc, const double boxsize,
		    const double k_min, const double k_max, const double dk,
		    const int mas_correction)
{
  const double fac= 2.0*M_PI/boxsize;
  //  const double sin_fac= 0.5*boxsize/nc;
  const int nbin= (int) round((k_max - k_min)/dk);

  init_wcorr(nc, mas_correction);
  
  //
  // Allocate memory
  //
  int* const nmodes= (int*) calloc(sizeof(int), nbin);
  double* const kmean= (double*) calloc(sizeof(double), nbin);
  double* const P= (double*) calloc(sizeof(double), nbin);
  //double dx= boxsize/nc;
  const size_t ncz= nc/2 + 1;

  for(int ix=0; ix<nc; ++ix) {
   double kx= ix <= nc/2 ? fac*ix : fac*(ix - nc);
   double w_x= wcorr[ix];

   for(int iy=0; iy<nc; ++iy) {
    double ky= iy <= nc/2 ? fac*iy : fac*(iy - nc);
    double w_xy= w_x*wcorr[iy];
       
    int iz0 = !(kx > 0.0 || (kx == 0.0 && ky > 0.0));

    for(int iz=iz0; iz<nc/2+1; ++iz) {
      double kz= fac*iz;
      double w_xyz= w_xy*wcorr[iz];
      
      size_t index= (ix*nc + iy)*ncz + iz;
      double d_re= fk[index][0];
      double d_im= fk[index][1];
      
      double k2= kx*kx + ky*ky + kz*kz;
      double k= sqrt(k2);

      int i= (int) floor((k - k_min)/dk);
      if(0 <= i && i < nbin) {
	nmodes[i]++;
	kmean[i] += k;
	P[i] += (d_re*d_re + d_im*d_im)/w_xyz;
      }
    }
   }
  }

  FILE* fp= fopen(filename, "w"); assert(fp);
  for(int i=0; i<nbin; ++i) {
    if(nmodes[i] > 0) {
      //double k= k_min + (i + 0.5)*dk;
      fprintf(fp, "%e %.15e %d\n",
	     kmean[i]/nmodes[i],
	     P[i]/nmodes[i],
	     nmodes[i]);
    }
  }
  fclose(fp);

  msg_printf(msg_info, "%s written.\n", filename);
}


