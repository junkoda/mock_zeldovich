#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>


static double Om= -1.0;
static double D0= 0.0;  // normalization factor

//static double Omega, OmegaLambda;
//static double growth_factor(const double a, const double a_normalization=1.0);

static double growth_int(double a, void* param);
static double growth(const double a);

void growth_init(const double omega_m)
{
  Om= omega_m;
  D0= growth(1.0);
}


double growth_D(const double a)
{
  // Growth factor D(a)
  assert(Om >= 0.0 && D0 > 0.0); // check growth_init() has done
  return growth(a)/D0;
}

double growth_f(const double a)
{
  assert(Om >= 0.0 && D0 > 0.0); // check growth_init() has done
  
  // Growth rate f=dlnD/dlna
  double hubble_a= sqrt(Om/(a*a*a) + 1.0 - Om);
  double d= growth(a)/D0;

  double f_ex= 1.0/(d*D0*a*a*hubble_a*hubble_a)
                 - 1.5*Om/(hubble_a*hubble_a*a*a*a);

  return f_ex;
}



double growth_int(double a, void* param)
{
  return pow(a/(Om + (1.0 - Om)*a*a*a), 1.5);
}

double growth(const double a)
{
  // Compute integral of D_tilda(a)
  // growth_factor = D_tilda(a)/D_tilda(1.0)
  
  const double hubble_a= sqrt(Om/(a*a*a) + 1.0 - Om);

  const int worksize= 1000;
  double result, abserr;

  gsl_integration_workspace *workspace=
    gsl_integration_workspace_alloc(worksize);

  gsl_function F;
  F.function = &growth_int;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, worksize, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}
