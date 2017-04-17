#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <boost/program_options.hpp>
//#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_integration.h>

#include "power.h"
#include "growth.h"
#include "corr.h"
#include "sincos_integ.h"

double zeldovich_power_realspace(const double k,
			       CorrelationFunction const * const psi0,
			       CorrelationFunction const * const psi2,
			       const int nmu);

//gsl_integration_cquad_workspace * w;
//static const int nwork= 100;
//static double sigma_v2;
//static const double dkr= 1.0;

using namespace std;
using namespace boost::program_options;

static void setup_r(const double rmax, const int n, vector<double>& v);
  
int main(int argc, char* argv[])
{
  //
  // command-line options (Boost program_options)
  //
  options_description opt("theory_zeldovich [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "filename")
    ("kmin", value<double>()->default_value(0.01), "k_min output")
    ("kmax", value<double>()->default_value(1.0), "k_max output")
    ("rmax", value<double>()->default_value(20000.0), "r maximum")
    ("nr", value<int>()->default_value(10000), "#bin in r integral")
    ("nmu", value<int>()->default_value(1000), "#bin in mu integral")
    ("a", value<double>()->default_value(1.0, "1"), "scale factor")
    ("omega_m", value<double>()->default_value(0.308, "0.308"), "omega_m")
    ("print-psi", "print correlation function Psi0(r) Psi1(r)")
    ("test", "test numerical accuracy of Fourier integral")
    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt; 
    return 0;
  }

  const string filename= vm["filename"].as<string>();
  const double rmax= vm["rmax"].as<double>();
  const int nr= vm["nr"].as<int>();
  const int nmu= vm["nmu"].as<int>();

  const double a= vm["a"].as<double>();
  const double omega_m= vm["omega_m"].as<double>();

  printf("# nrbin= %d, nmubin= %d, rmax= %e\n", nr, nmu, rmax);

  // Setup redshift growth
  growth_init(omega_m);
  const double f= growth_f(a); 
  const double growth= growth_D(a);
  printf("# D= %f, f= %f\n", growth, f);
    

  // Setup power spectrum
  PowerSpectrum* const ps= power_alloc(filename.c_str(), 0);
  power_rescale(ps, growth);

  // Setup correlation_functions
  vector<double> vr;
  setup_r(rmax, nr, vr);

  vector<CorrelationFunction*> psi
    = corr_compute_psi(ps, vr);
  CorrelationFunction const * const psi0= psi[0];
  CorrelationFunction const * const psi2= psi[1];
  

  if(vm.count("print-psi")) {
    corr_print(psi0, psi2);
    return 0;
  }
  else if(vm.count("test")) {
    corr_fourier_test(psi0, psi2);
    return 0;
  }

  const double k_min= vm["kmin"].as<double>();
  const double k_max= vm["kmax"].as<double>();
  
  for(int i=0; i<ps->n; ++i) {
    double k= ps->k[i];
    if(k_min <= k && k <= k_max) {
      double Pz= zeldovich_power_realspace(k, psi0, psi2, nmu);
      printf("%e %e %e\n", k, ps->P[i], Pz);
    }
  }

  // Column 1: k
  // Column 2: P_linear
  // Column 3: P_Zeldovich

  
  return 0;
}

void setup_r(const double rmax, const int n, vector<double>& v)
{
  const double logrmin= log(0.001);
  const double logrmax= log(rmax);
  for(int i=0; i<n; ++i) {
    double r= exp(logrmin + (logrmax - logrmin)/(n-1)*i);
    v.push_back(r);
  }
}


double zeldovich_power_realspace(const double k,
				 CorrelationFunction const * const psi0,
				 CorrelationFunction const * const psi2,
				 const int nmu)
{
  //
  // 2D Fourier transform Psi back to P(k) in linear limit
  //
  // P(k) = 4pi k^2 \int_0^1 dcosO \int dr r^2 dr
  //          [ 1/3 \Psi_0(r) - (cos^2 O - 1/3) \Psi_2(r) ]
  //
  const size_t n= psi0->n; assert(n > 0);
  const double sigma_v2 = psi0->xi[0]/3.0;

  const double dmu= 1.0/nmu;

  //for(double k= 0.01; k<1.0; k+=0.01) {
  double integ= 0.0;
  
  for(int imu=0; imu<nmu; ++imu) {  // dcosO integral
    double mu= (imu + 0.5)/nmu;     // mu = cosO
    double cos2 = mu*mu - 1.0/3.0;  // cos^2 O - 1/3
    
    double integ_krc= 0.0;
    for(size_t i=1; i<n; ++i) {
      // krc = k r cosO
      // r^2 dr = 1/(k cos)^3 O krc^2 dkr
      double krc1= k*psi0->r[i-1]*mu;
      double krc2= k*psi0->r[i]*mu;
      
      double sinkrc1= sin(krc1);
      double coskrc1= cos(krc1);
      
      double sinkrc2= sin(krc2);
      double coskrc2= cos(krc2);
      
      double r1 = psi0->r[i-1];
      double r2 = psi0->r[i];
      
      double xi1= exp(k*k*((psi0->xi[i-1] - psi0->xi[0])/3.0
			   - cos2*psi2->xi[i-1]))
	            - exp(-k*k*psi0->xi[0]/3.0);
      double xi2= exp(k*k*((psi0->xi[i  ] - psi0->xi[0])/3.0
				 - cos2*psi2->xi[i]))
	            - exp(-k*k*psi0->xi[0]/3.0);

      // a0 + a1*kr
      double a1= (xi2 - xi1)/(krc2 - krc1);
      double a0= (xi1*krc2 - xi2*krc1)/(krc2 - krc1);
      
      integ_krc +=  a1*(cos_integ3(krc2, sinkrc2, coskrc2) -
			cos_integ3(krc1, sinkrc1, coskrc1)) +
	          + a0*(cos_integ2(krc2, sinkrc2, coskrc2) -
			cos_integ2(krc1, sinkrc1, coskrc1));
    }
    integ += integ_krc*dmu/(mu*mu*mu);
  }    
  return 4.0*M_PI/(k*k*k)*integ;
}
