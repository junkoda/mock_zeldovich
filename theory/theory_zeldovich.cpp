#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <boost/program_options.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

#include "power.h"
#include "growth.h"
#include "corr.h"
#include "sincos_integ.h"

gsl_integration_cquad_workspace * w;
static const int nwork= 100;
static double sigma_v2;
static const double dkr= 1.0;

//#include "power_spectrum.h"
//#include "correlation_functions.h"
//

double integrate_1(const double k, const double r[], const int n);
double integrate_dd(const double k, const double r[], const int n);
double integrate_du(const double k, const double r[], const int n);
double integrate_du2(const double k, const double r[], const int n);
  
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
    ("rmax", value<double>()->default_value(1000.0), "r maximum")
    ("nbin", value<int>()->default_value(1000), "nbin")
    ("z", value<double>()->default_value(0.0, "0"), "redshift")
    ("omega_m", value<double>()->default_value(0.308, "0.308"), "omega_m")
    ("print-psi", "print correlation function Psi0(r) Psi1(r)")
    ("test1", "test numerical accuracy of 1D Fourier integral")
    ("test2", "test numerical accuracy of 2D Fourier integral")
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
  const int nbin= vm["nbin"].as<int>();

  const double z= vm["z"].as<double>();
  const double a= 1.0/(1+z);
  const double omega_m= vm["omega_m"].as<double>();

  printf("# nrbin= %d, rmax= %e\n", nbin, rmax);

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
  setup_r(rmax, nbin, vr);
  
  corr_init(&vr.front(), vr.size(), ps);

  corr_integrate_psi_hybrid();

  if(vm.count("print-psi")) {
    corr_print();
    return 0;
  }
  else if(vm.count("test1")) {
    corr_fourier_test1();
    return 0;
  }
  else if(vm.count("test2")) {
    corr_fourier_test2();
    return 0;
  }



  /*
  double sigma_v= f*power_sigmav(ps);
  sigma_v2= sigma_v*sigma_v;
  printf("# sigma_v= %f\n", f*sigma_v);
  
  w= gsl_integration_cquad_workspace_alloc(nwork);

  cerr << "start Pi integration\n";
  for(double k= 0.01; k<1.0; k+=0.01) {
    double Pi_dd= integrate_dd(k, &vr.front(), vr.size());
    double Pi_du= integrate_du(k, &vr.front(), vr.size());
    double Pi_du2= integrate_du2(k, &vr.front(), vr.size());

    // debug!!! integrate_1 not working
    //double Pi_1= integrate_1(k, &vr.front(), vr.size());

    printf("%e %e %e %e %e\n", k, Pi_dd, Pi_du, Pi_du2, 0.0);
    //printf("%e %e %e %e %e\n", k, 0.0, 0.0, 0.0, Pi_1);
  }

  gsl_integration_cquad_workspace_free(w);
  */
  
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

/*
double integrand_cos(double mu, void* vparams)
{
  double* params = (double *) vparams;
  double k= params[0];
  double kr= params[1];
  double xi_uu0= params[2];
  double xi_uu2= params[3];
  double P2= 1.5*mu*mu - 0.5;

  return cos(kr*mu)*exp(-k*k*(sigma_v2 - xi_uu0 - xi_uu2*P2));
  //return cos(kr*mu)*exp(-k*k*ff*sigma_v2);
  //return cos(kr*mu); // debug!!
}

double integrand_cos1(double mu, void* vparams)
{
  double* params = (double *) vparams;
  double k= params[0];
  double kr= params[1];
  double xi_uu0= params[2];
  double xi_uu2= params[3];
  double P2= 1.5*mu*mu - 0.5;

  return cos(kr*mu)*(k*k*(xi_uu0 + xi_uu2*P2));
  //return cos(kr*mu)*(exp(-k*k*(sigma_v2 - xi_uu0 - xi_uu2*P2)) - 1.0);
  //return cos(kr*mu)*exp(-k*k*ff*sigma_v2);
  //return cos(kr*mu); // debug!!
}

double integrand_mu_sin(double mu, void* vparams)
{
  double* params = (double *) vparams;
  double k= params[0];
  double kr= params[1];
  double xi_uu0= params[2];
  double xi_uu2= params[3];
  double P2= 1.5*mu*mu - 0.5;
  
  return mu*sin(kr*mu)*exp(-k*k*(sigma_v2 - xi_uu0 - xi_uu2*P2));
  //return mu*sin(kr*mu); // debug!
}

double integrand_mu2_cos(double mu, void* vparams)
{
  double* params = (double *) vparams;
  double k= params[0];
  double kr= params[1];
  double xi_uu0= params[2];
  double xi_uu2= params[3];
  double P2= 1.5*mu*mu - 0.5;
  
  return mu*mu*cos(kr*mu)*exp(-k*k*(sigma_v2 - xi_uu0 - xi_uu2*P2));
  //return mu*mu*cos(kr*mu); // debug!
}

double integrate_dd(const double k, const double r[], const int n)
{  
  gsl_function F;
  double params[4];
  
  F.function = &integrand_cos;
  F.params   = params;

  double result, error;
  double r0= r[0];
  double f0= 0.0;
  
  double integ= 0.0;
    
  for(int i=0; i<n; ++i) {
    double r1= r[i];
    params[0]= k;
    params[1]= k*r1;
    params[2]= corr_uu0(i);
    params[3]= corr_uu2(i);

    gsl_integration_cquad(&F, 0.0, 1.0, 0.0, 1.0e-5, w, &result, NULL, NULL);

    double f1= result*corr_dd(i);

    integ += 0.5*(r0*r0*f0 + r1*r1*f1)*(r1 - r0);

    r0= r1;
    f0= f1;
  }

  const double fac= 4.0*M_PI;
  return fac*integ;
}

double integrate_du(const double k, const double r[], const int n)
{  
  gsl_function F;
  double params[4];
  
  F.function = &integrand_mu_sin;
  F.params   = params;

  double result, error;
  double r0= r[0];
  double f0= 0.0;
  
  double integ= 0.0;
    
  for(int i=0; i<n; ++i) {
    double r1= r[i];
    params[0]= k;
    params[1]= k*r1;
    params[2]= corr_uu0(i);
    params[3]= corr_uu2(i);

    gsl_integration_cquad(&F, 0.0, 1.0, 0.0, 1.0e-5, w, &result, NULL, NULL);

    double f1= result*corr_du1(i);

    integ += 0.5*(r0*r0*f0 + r1*r1*f1)*(r1 - r0);

    r0= r1;
    f0= f1;
  }

  const double fac= 4.0*M_PI;
  return fac*integ*2.0*k;
}

double integrate_du2(const double k, const double r[], const int n)
{  
  gsl_function F;
  double params[4];
  
  F.function = &integrand_mu2_cos;
  F.params   = params;

  double result, error;
  double r0= r[0];
  double f0= 0.0;
  
  double integ= 0.0;
    
  for(int i=0; i<n; ++i) {
    double r1= r[i];
    params[0]= k;
    params[1]= k*r1;
    params[2]= corr_uu0(i);
    params[3]= corr_uu2(i);

    gsl_integration_cquad(&F, 0.0, 1.0, 0.0, 1.0e-5, w, &result, NULL, NULL);

    double xi_du1= corr_du1(i);
    double f1= result*xi_du1*xi_du1;

    integ += 0.5*(r0*r0*f0 + r1*r1*f1)*(r1 - r0);

    r0= r1;
    f0= f1;
  }

  const double fac= 4.0*M_PI;
  return fac*integ*k*k;
}


double integrate_1(const double k, const double r[], const int n)
{  
  gsl_function F;
  double params[4];
  
  F.function = &integrand_cos1;
  F.params   = params;

  double result, error;
  double r0= r[0];
  double f0= 0.0;
  
  double integ= 0.0;
    
  for(int i=0; i<n; ++i) {
    double r1= r[i];
    params[0]= k;
    params[1]= k*r1;
    params[2]= corr_uu0(i);
    params[3]= corr_uu2(i);

    gsl_integration_cquad(&F, 0.0, 1.0, 0.0, 1.0e-5, w, &result, NULL, NULL);

    double f1= result;

    integ += 0.5*(r0*r0*f0 + r1*r1*f1)*(r1 - r0);

    r0= r1;
    f0= f1;
  }

  const double fac= 4.0*M_PI;
  return fac*integ*k*k;
}
*/
