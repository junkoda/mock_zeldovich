//
// Analytical functions in cosmology, e.g. linear growth rate
// Flat Lambda CDM assumed
//

#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "msg.h"
#include "const.h"
#include "cosmology.h"
#include "error.h"

namespace {
  double omega_m0;
  double growth_normalisation;

  void check_initialisation();
  double growth_integrand(double a, void* param);
  double growth_unnormalised(const double a);
}

void cosmology_init(const double omega_m0_)
{
  omega_m0= omega_m0_;

  msg_printf(msg_info, "Cosmology initialised with omega_m = %.7f\n", omega_m0);

  growth_normalisation= 1.0/growth_unnormalised(1.0); // D_growth=1 at a=1
}

double cosmology_D_growth(const double a)
{
  check_initialisation();
  // Linear growth factor D
  if(a == 0.0) return 0.0;
  
  return growth_normalisation*growth_unnormalised(a);
}

double cosmology_D2_growth(const double a, const double D)
{
  // 2nd-order growth factor D2
  if(a == 0.0) return 0.0;
  
  return -3.0/7.0*D*D*pow(cosmology_omega(a), -1.0/143.0);
}

double cosmology_Dv_growth(const double a, const double D)
{
  if(a == 0.0) return 1.0;

  double H= cosmology_hubble_function(a);
  double f= cosmology_f_growth_rate(a);
  
  return a*a*D*H*f;
}

double cosmology_D2v_growth(const double a, const double D2)
{
  double H= cosmology_hubble_function(a);
  double f= cosmology_f_growth_rate(a);

  return 2.0*a*a*D2*H*f;
}

double cosmology_D2a_growth(const double D1, const double D2)
{
  return D2 - D1*D1;
}

double cosmology_f_growth_rate(const double a)
{
  check_initialisation();

  if(a == 0.0) return 1.0;
  
  // Linear growth rate f=dlnD/dlna
  const double d_un= growth_unnormalised(a);
  const double hf= cosmology_hubble_function(a);

  return 1.0/(d_un*a*a*hf*hf) - 1.5*omega_m0/(hf*hf*a*a*a);   
}  


void cosmology_growth(const double a,
		      double* const D_result, double* const f_result)
{
  check_initialisation();
  // Both linear growth factor D(a) and growth rate f=dlnD/dlna

  if(a == 0.0) {
    *D_result= 0.0;
    *f_result= 1.0;
    return;
  }
  
  const double d_un= growth_unnormalised(a);
  const double hf= cosmology_hubble_function(a);

  *D_result= growth_normalisation*d_un;
  *f_result= 1.0/(d_un*a*a*hf*hf) - 1.5*omega_m0/(hf*hf*a*a*a);   
}

double cosmology_hubble_function(const double a)
{
  // H/H0= sqrt(Omega_m0*a^-3 + Omega_Lambda)
  return sqrt(omega_m0/(a*a*a) + (1 - omega_m0));
}

double cosmology_omega(const double a)
{
  // Omega_m(a)
  check_initialisation();
  return omega_m0/(omega_m0 + (1 - omega_m0)*(a*a*a));
}

double cosmology_omega_m()
{
  check_initialisation();
  return omega_m0;
}

double cosmology_rho_m()
{
  check_initialisation();
  return omega_m0*c::rho_crit_0;
}


namespace {

void check_initialisation()
{
  // Check if this module is initilised
  if(growth_normalisation == 0.0) {
    msg_printf(msg_fatal, "Error: cosmology module not initialised.\n");
    throw RuntimeError();
  }
}
    

double growth_integrand(double a, void* param)
{
  // sqrt[(a*H/H0)^3]
  //return pow(a/(omega_m0 + (1 - Omega)*a*a*a), 1.5);
  const double aHinv= 1.0/sqrt(omega_m0/a + (1 - omega_m0)*(a*a));
  return aHinv*aHinv*aHinv;
}


double growth_unnormalised(const double a)
{
  // D(a) \propto \int_0^a (a H(a)/H0)^-3 da
  const size_t worksize= 1000;

  gsl_integration_workspace *workspace=
    gsl_integration_workspace_alloc(worksize);

  gsl_function F;
  F.function = &growth_integrand;

  double result, abserr;
  gsl_integration_qag(&F, 0, a, 0, 0.5e-8, worksize, GSL_INTEG_GAUSS41,
		      workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return cosmology_hubble_function(a) * result;
}

}
