#ifndef SINCOS_INTEG_H
#define SINCOS_INTEG_H 1

static inline double sin_integ0(const double kr, const double sinkr, const double coskr)
{
  // \int sin(kr) dkr = k0sin_integ(kr)
  return -coskr;
}

static inline double sin_integ1(const double kr, const double sinkr, const double coskr)
{
  // \int kr sin(kr) dkr = k1sin_integ(kr)
  return sinkr - kr*coskr;
}
static inline double sin_integ2(const double kr, const double sinkr, const double coskr)
{
  // \int kr^2 sin(kr) dkr = k2sin_integ(kr)
  return (2.0 - kr*kr)*coskr + 2.0*kr*sinkr;
}

static inline double cos_integ0(const double x, const double sinx, const double cosx)
{
  // \int cos(x) dx
  return sinx;
}

static inline double cos_integ1(const double x, const double sinx, const double cosx)
{
  // \int x cos(x) dx
  return x*sinx + cosx;
}

static inline double cos_integ2(const double x, const double sinx, const double cosx)
{
  // \int x^2 cos(x) dx
  return (x*x - 2.0)*sinx + 2.0*x*cosx;
}

static inline double cos_integ3(const double x, const double sinx, const double cosx)
{
  // \int x^3 cos(x) dx
  return (3.0*x*x - 6.0)*cosx + (x*x*x - 6.0*x)*sinx;
}





#endif
