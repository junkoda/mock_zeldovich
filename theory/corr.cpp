//
// Compute correlation function <Psi_i(p) Psi_j(q)>
// of Zeldovich displacement Psi_i(p) from power spectrum P(k)
//
// <Psi_i(p) Psi_j(q)>
// = Psi_parallel(r) r^_i r^_j + Psi_perp(r) [ delta_ij - r^_i r^_j ]
// = [ Psi0(r) - 2 Psi1(r) ] r^_i r^_j + Psi_1 [ delta_ij - r^_i r^_j ]
//
// intergrate_psi_hybrid() computes psi0 and psi1

#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include "corr.h"
#include "power.h"
#include "sincos_integ.h"

using namespace std;

CorrelationFunction* corr_alloc(vector<double>& r);

//
// Global variables
//
static PowerSpectrum const *p;
//static CorrelationFunction *psi0, *psi2;
//static gsl_integration_workspace *w;

//
// Static functions
//
void integrate_psi_trapezoidal_r(const double r, double * const psi);
void integrate_psi_hybrid_r(const double r, double * const psi);


vector<CorrelationFunction*>
corr_compute_psi(PowerSpectrum const * const ps, vector<double>& vr)
{
  p= ps;

  CorrelationFunction* psi0= corr_alloc(vr);
  CorrelationFunction* psi2= corr_alloc(vr);


  
  double const * const r= psi0->r;
  const size_t n= psi0->n;

  for(size_t i=0; i<n; ++i) {
    double psi[3];

    integrate_psi_hybrid_r(r[i], psi);
    //integrate_psi_trapezoidal_r(r[i], psi); // DEBUG!!

    psi0->xi[i]= psi[0];
    psi2->xi[i]= psi[2];
  }

  vector<CorrelationFunction*> psi;
  psi.push_back(psi0);
  psi.push_back(psi2);

  return psi;
}

CorrelationFunction* corr_alloc(vector<double>& vr)
{
  size_t n= vr.size();
  
  CorrelationFunction* corr=
    (CorrelationFunction*) malloc(sizeof(CorrelationFunction));

  
  corr->n= n;
  corr->r= (double*) malloc(sizeof(double)*n);
  corr->xi= (double*) malloc(sizeof(double)*n);

  for(int i=0; i<n; ++i)
    corr->r[i]= vr[i];

  return corr;
}

void corr_print(CorrelationFunction const * const psi0,
		CorrelationFunction const * const psi2)
{
  printf("# r psi0 psi2\n");

  const double mu = 1.0;
  const size_t n= psi0->n;
  for(size_t i=0; i<n; ++i) {
    double e1 = 1.0/3.0*(psi0->xi[i] - psi0->xi[0]);
    double e2 = (3.0*mu*mu - 1.0)*psi2->xi[i];
    
    printf("%le %le %le %le %le %le\n",
	   psi0->r[i], psi0->xi[i], psi2->xi[i],
	   e1, e2, exp(e1 + e2));
  }
  // Column 1: r
  // Column 2: Psi0(r)
  // COlumn 3: Psi2(r)
}

void integrate_psi_trapezoidal_r(const double r, double * const psi)
{
  const int n= p->n;
  const double ktr= 0.005;

  double xi0= 0.0, xi2= 0.0;
  for(int i=1; i<n; ++i) {
    //Multiply Ptt(k) by k^2 for k<ktr to avoid divergent psi(r)*r^2
    double k1= p->k[i-1];
    double k2= p->k[i];

    double fac1= k1 < ktr ? (k1/ktr)*(k1/ktr) : 1.0;
    double fac2= k2 < ktr ? (k2/ktr)*(k2/ktr) : 1.0;

    double P1= fac1*p->P[i-1];
    double P2= fac2*p->P[i];

    double k1r= k1*r;
    double k2r= k2*r;

    // += P(k) j0(kr) dkr

    xi0 += 0.5*(gsl_sf_bessel_j0(k1r)*P1 +
		gsl_sf_bessel_j0(k2r)*P2)*(k2r - k1r);

    // += P(k) j2(kr) dkr
    xi2 += 0.5*(gsl_sf_bessel_j2(k1r)*P1 +
		gsl_sf_bessel_j2(k2r)*P2)*(k2r - k1r);
  }

  const double fac= 1.0/(2.0*M_PI*M_PI*r);

  psi[0]= fac*xi0;
  psi[2]= fac*xi2;
}

  
  
void integrate_psi_hybrid_r(const double r, double * const psi)
{
  // Prerequisit:
  //     PowerSpectrum p
  // Input:
  //     r [1/h Mpc]
  // Ouput:
  //     psi[0] = 1/(2 pi^2) \int dk P(k) j0(kr)
  //     psi[1] = 1/(2 pi^2) \int dk P(k) j1(kr)/(kr)
  //     psi[2] = 1/(2 pi^2) \int dk P(k) j2(kr)
  //
  // where j0, j1 are spherical Bessel functions:
  //     j0(kr) = sin(kr)/kr
  //     j1(kr) = sin(kr)/(kr)^2 - cos(kr)/(kr)
  
  double xi0= 0.0, xi2= 0.0;

  const int n= p->n;
  const double ktr= 0.005;
    
  for(int i=1; i<n; ++i) {
    //Multiply Ptt(k) by k^2 for k<ktr to avoid divergent psi(r)*r^2
    double k1= p->k[i-1];
    double k2= p->k[i];

    double fac1= k1 < ktr ? (k1/ktr)*(k1/ktr) : 1.0;
    double fac2= k2 < ktr ? (k2/ktr)*(k2/ktr) : 1.0;

    double P1= fac1*p->P[i-1];
    double P2= fac2*p->P[i];

    double k1r= k1*r;
    double k2r= k2*r;

    double sink1r= sin(k1r);
    double cosk1r= cos(k1r);

    double sink2r= sin(k2r);
    double cosk2r= cos(k2r);


    // \int P(k) j0(kr) dkr

    if(k1r < 1.0) {
      // Trapezoidal for smooth function
      xi0 += 0.5*(gsl_sf_bessel_j0(k1r)*P1 +
		  gsl_sf_bessel_j0(k2r)*P2)*(k2r - k1r);
    }
    else {
      // Integration of oscilating function
      // Linear interpolation P(k)/(kr) = a1_inv1*kr + a0_inv1

      double a1_inv1= (P2/k2r - P1/k1r)/(k2r - k1r);
      double a0_inv1= (P1/k1r*k2r - P2/k2r*k1r)/(k2r - k1r);

      xi0 += a0_inv1*(sin_integ0(k2r, sink2r, cosk2r) -
		      sin_integ0(k1r, sink1r, cosk1r))
	   + a1_inv1*(sin_integ1(k2r, sink2r, cosk2r) -
		      sin_integ1(k1r, sink1r, cosk1r));

    }

    // \int P(k) j2(kr) dkr
    if(k1r < 1.0) {
      // Trapezoidal
      xi2 += 0.5*(gsl_sf_bessel_j2(k1r)*P1 +
		  gsl_sf_bessel_j2(k2r)*P2)*(k2r - k1r);
    }
    else {
      // analytic integration of (a1*kr + a0)*cos(kr) or sin(kr)
      // Linear interpolation P(k)/(kr)^2 = a1_inv1*kr + a0_inv1
      double a1_inv1= (P2/k2r - P1/k1r)/(k2r - k1r);
      double a0_inv1= (P1/k1r*k2r - P2/k2r*k1r)/(k2r - k1r);

      // - sin(kr)/(kr) P(k) dkr
      xi2 -= a0_inv1*(sin_integ0(k2r, sink2r, cosk2r) -
		      sin_integ0(k1r, sink1r, cosk1r))
	   + a1_inv1*(sin_integ1(k2r, sink2r, cosk2r) -
		      sin_integ1(k1r, sink1r, cosk1r));

      // Linear interpolation P(k)/(kr)^2 = a1_inv2*kr + a0_inv2
      double k1r2= k1r*k1r;
      double k2r2= k2r*k2r;

      double a1_inv2= (P2/k2r2 - P1/k1r2)/(k2r - k1r);
      double a0_inv2= (P1/k1r2*k2r - P2/k2r2*k1r)/(k2r - k1r);

      // - 3 \int cos(kr)/(kr)^2 dkr
      xi2 -= 3.0*(a0_inv2*(  cos_integ0(k2r, sink2r, cosk2r)
			   - cos_integ0(k1r, sink1r, cosk1r))
		+ a1_inv2*(  cos_integ1(k2r, sink2r, cosk2r)
			   - cos_integ1(k1r, sink1r, cosk1r)));


      // Linear interpolation P(k)/(kr)^3 = a1_inv3*kr + a0_inv3
      double k1r3= k1r2*k1r;
      double k2r3= k2r2*k2r;

      double a1_inv3= (P2/k2r3 - P1/k1r3)/(k2r - k1r);
      double a0_inv3= (P1/k1r3*k2r - P2/k2r3*k1r)/(k2r - k1r);

      // 3 \int sin(kr)/(kr)^3 dkr
      xi2 += 3.0*(a0_inv3*(  sin_integ0(k2r, sink2r, cosk2r)
			   - sin_integ0(k1r, sink1r, cosk1r))
		+ a1_inv3*(  sin_integ1(k2r, sink2r, cosk2r)
			   - sin_integ1(k1r, sink1r, cosk1r)));


    }
  }

  const double fac= 1.0/(2.0*M_PI*M_PI*r);

  psi[0]= fac*xi0;
  psi[2]= fac*xi2;
}

void corr_fourier_test(CorrelationFunction const * const psi0,
		       CorrelationFunction const * const psi2)
{  
  // Fourier transform psi(r) back to P(k)
  printf("# fourier_test\n");
  const size_t n= psi0->n; assert(n > 0);

  for(double k= 0.01; k<1.0; k+=0.01) {
    double integ0= 0.0;
    double integ2= 0.0;

    for(size_t i=1; i<n; ++i) {
      double kr1= k*psi0->r[i-1];
      double kr2= k*psi0->r[i];

      double sinkr1= sin(kr1);
      double coskr1= cos(kr1);

      double sinkr2= sin(kr2);
      double coskr2= cos(kr2);

      //
      // 4pi \int psi0(r) j0(kr) r^2 dr = P(k)/k^2
      //
      double xi1= psi0->xi[i-1];
      double xi2= psi0->xi[i];

      // psi0(r) ~ a0 + a1*kr
      double a1= (xi2 - xi1)/(kr2 - kr1);
      double a0= (xi1*kr2 - xi2*kr1)/(kr2 - kr1);

      // integrate sin(kr)*(kr)*(a0 + a1*kr)
      integ0 +=  a1*(sin_integ2(kr2, sinkr2, coskr2) -
		     sin_integ2(kr1, sinkr1, coskr1)) +
	       + a0*(sin_integ1(kr2, sinkr2, coskr2) -
		     sin_integ1(kr1, sinkr1, coskr1));

      //
      // 4pi \int psi1(r) j2(kr) r^2 dr = P(k)/k^2
      //
      // (kr)^2 j2(kr) = [3.0/(kr) - (kr)] sin(kr) - 3 cos(kr)

      // psi2(r) ~ a0 + a1*kr
      xi1= psi2->xi[i-1];
      xi2= psi2->xi[i];

      a1= (xi2 - xi1)/(kr2 - kr1);
      a0= (xi1*kr2 - xi2*kr1)/(kr2 - kr1);

      // - 3 cos(kr) [a0 + a1*kr] dk
      integ2 -= 3.0*(a0*(cos_integ0(kr2, sinkr2, coskr2) -
	                 cos_integ0(kr1, sinkr1, coskr1))
	           + a1*(cos_integ1(kr2, sinkr2, coskr2) -
	                 cos_integ1(kr1, sinkr1, coskr1)));

      // - sin(kr) (kr) [a0 + a1*kr] dk
      integ2 -= a0*(sin_integ1(kr2, sinkr2, coskr2) -
	            sin_integ1(kr1, sinkr1, coskr1))
	      + a1*(sin_integ2(kr2, sinkr2, coskr2) -
	            sin_integ2(kr1, sinkr1, coskr1));


      // psi2(r)/kr ~ a0 + a1*kr
      a1= (xi2/kr2 - xi1/kr1)/(kr2 - kr1);
      a0= (xi1/kr1*kr2 - xi2/kr2*kr1)/(kr2 - kr1);

      // 3/kr^3 sin(kr) dkr ~ -3 sin(kr) [a0 + a1*kr] dk
      integ2 += 3.0*(a0*(sin_integ0(kr2, sinkr2, coskr2) -
	                 sin_integ0(kr1, sinkr1, coskr1))
		   + a1*(sin_integ1(kr2, sinkr2, coskr2) -
			 sin_integ1(kr1, sinkr1, coskr1)));
      

    }

    const double fac= 4.0*M_PI/k;

    printf("%e %e %e\n",
	   k,
	   fac*integ0,
	   fac*integ2);
    // fourier transform back
    // Column 1: k
    // Column 2: P(k) transformed back from psi0
    // Column 3: P(k) transformed back from psi2
  }
}


/*
double corr_dd(const int i)
{
  return xi_dd->xi[i];
}

double corr_du1(const int i)
{
  return xi_du1->xi[i];
}

double corr_uu0(const int i)
{
  return xi_uu0->xi[i];
}

double corr_uu2(const int i)
{
  return xi_uu2->xi[i];
}

double corr_dd(const double r)
{
  if(r < xi_dd->r[0] || r > xi_dd->r[xi_dd->n-1]) return 0.0;
  return gsl_interp_eval(xi_dd->interp, xi_dd->r, xi_dd->xi, r, xi_dd->acc);
}

double corr_du1(const double r)
{
  if(r < xi_du1->r[0] || r > xi_du1->r[xi_du1->n-1]) return 0.0;
  return gsl_interp_eval(xi_du1->interp, xi_du1->r, xi_du1->xi, r, xi_du1->acc);
}

double corr_uu0(const double r)
{
  if(r < xi_uu0->r[0] || r > xi_uu0->r[xi_uu0->n-1]) return 0.0;  
  return gsl_interp_eval(xi_uu0->interp, xi_uu0->r, xi_uu0->xi, r, xi_uu0->acc);
}

double corr_uu2(const double r)
{
  if(r < xi_uu2->r[0] || r > xi_uu2->r[xi_uu2->n-1]) return 0.0;
  return gsl_interp_eval(xi_uu2->interp, xi_uu2->r, xi_uu2->xi, r, xi_uu2->acc);
}

CorrelationFunction const * get_xi_dd()
{
  return xi_dd;
}

CorrelationFunction const * get_xi_du1()
{
  return xi_du1;
}

*/

