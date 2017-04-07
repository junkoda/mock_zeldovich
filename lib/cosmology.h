#ifndef COSMOLOGY_H
#define COSMOLOGY_H 1

void   cosmology_init(const double omega_m0);

double cosmology_D_growth(const double a);
double cosmology_D2_growth(const double a, const double D);
double cosmology_Dv_growth(const double a, const double D);
double cosmology_D2v_growth(const double a, const double D2);
double cosmology_D2a_growth(const double D1, const double D2);

void   cosmology_growth(const double a, double* const D, double* const f);

double cosmology_f_growth_rate(const double a);

double cosmology_hubble_function(const double a);
double cosmology_omega(const double a);

double cosmology_omega_m();
double cosmology_rho_m();
#endif
