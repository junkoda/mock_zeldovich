//
// Definition of particle & particles data structure
//

#ifndef PARTICLE_H
#define PARTICLE_H 1

#include <stdint.h>
#include <vector>
#include "config.h"

struct Particle {
  Float  x[3];
  Float  v[3];
  Float  dx1[3];  // 1LPT (ZA) displacement
  Float  dx2[3];  // 2LPT displacement
  uint64_t id;      // Particle index 1,2,3...
  //uint64_t igrp;      // Group index
};

struct Pos {
  Float x[3];
};

class Particles {
 public:
  Particles(const int nc, const double boxsize);
  ~Particles();
  void update_np_total();
  
  Particle* p;
  std::vector<Particle>* pv;
  double a_x, a_v;
  Float3* force;

  size_t np_local, np_allocated;
  uint64_t np_total;
  double boxsize;
};

void particles_update_np_total(Particles* const particles);

#endif
