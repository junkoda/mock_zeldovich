#ifndef UTIL_H
#define UTIL_H 1

#include <string.h>
#include <assert.h>
#include "particle.h"

static inline size_t mbytes(size_t bytes)
{
  return bytes/(1024*1024);
}


static inline void periodic_wrapup_p(Particle& p, const Float boxsize)
{
  for(int k=0; k<3; k++) {
    if(p.x[k] < 0) p.x[k] += boxsize;
    if(p.x[k] >= boxsize) p.x[k] -= boxsize;

#ifdef CHECK
    assert(0 <= p.x[k] && p.x[k] < boxsize);
#endif
  }
}

static inline Float periodic_wrapup_x(Float x, const Float boxsize)
{
  if(x < 0) x += boxsize;
  if(x >= boxsize) x -= boxsize;

#ifdef CHECK
  assert(0 <= x && x < boxsize);
#endif

  return x;
}

bool util_stat(const char filename);

void util_periodic_wrapup(Particles* const particles);

#endif
