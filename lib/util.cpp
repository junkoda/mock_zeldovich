#include "util.h"

void util_periodic_wrapup(Particles* const particles)
{
  const size_t n= particles->np_local;
  Particle* const p= particles->p;
  const Float boxsize= particles->boxsize;
  
  for(size_t i=0; i<n; ++i)
    periodic_wrapup_p(p[i], boxsize);
}

