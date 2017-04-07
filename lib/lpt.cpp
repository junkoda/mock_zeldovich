//
// Lagrangian Perturbation Theory
//
// Generates random Gaussian pertubation using 2LPT
//

#include <iostream>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include "comm.h"
#include "msg.h"
#include "mem.h"
#include "config.h"
#include "cosmology.h"
#include "power.h"
#include "particle.h"
#include "fft.h"
#include "lpt.h"

using namespace std;

namespace {
  unsigned int* seedtable= 0;
  double boxsize;

  size_t nc= 0, ncz, nckz;
  size_t local_nx;
  size_t local_ix0;
  Float offset= 0.5f;

  FFT* fft_psi[3];    // Zeldovichi displacement Psi_i

  void set_seedtable(const int nc, gsl_rng* random_generator,
		     unsigned int* const stable);
}

static inline size_t i_x(int ix, int iy, int iz) {
  return (ix*nc + iy)*ncz + iz;
}

/*
static inline size_t i_k(int ikx, int iky, int ikz) {
  return (ikx*nc + iky)*nckz + ikz;
}
*/


void lpt_init(const int nc_, const double boxsize_, Mem* mem)
{
  // nc_: number of particles per dimension
  // boxsize_: box length on a side
  // Mem: Memory object for LPT, can be 0.
  //      If mem = 0, memory is allocated exclusively for LPT
  boxsize= boxsize_;
  nc= nc_;
  nckz= nc/2 + 1;

  msg_printf(msg_debug, "lpt_init(nc= %d, boxsize= %.1lf)\n", nc, boxsize);

  if(mem)
    mem->use_from_zero(0);
  
  for(int i=0; i<3; i++)
    fft_psi[i]= new FFT("Psi_i", nc, mem, 0);

  seedtable = (unsigned int *) malloc(nc*nc*sizeof(unsigned int)); assert(seedtable);

  // checks
  local_nx= fft_psi[0]->local_nx;
  local_ix0= fft_psi[0]->local_ix0;
  
  for(int i=0; i<3; i++) {
    assert(fft_psi[i]->nc == nc);
    assert(fft_psi[i]->local_nx == local_nx);
    assert(fft_psi[i]->local_ix0 == local_ix0);
  }
}

void lpt_free()
{
  for(int i=0; i<3; i++)
    delete fft_psi[i];

  free(seedtable);
  seedtable= 0;

  nc= 0;
}


void lpt_generate_psi_x(const unsigned long seed,
			double const * const Pk,
			const bool fix_amplitude)
{
  // Generates 1LPT (Zeldovich) displacements, Psi_k
  // from N-GenIC by Volker Springel

  // Pk is a 3D grid of power spectrum
  //
  msg_printf(msg_verbose, "Generating delta_k...\n");
  msg_printf(msg_info, "Random Seed = %lu\n", seed);

  for(int i=0; i<3; i++) assert(fft_psi[i]);
  
  complex_t* psi_k[]= {fft_psi[0]->fk, fft_psi[1]->fk, fft_psi[2]->fk};

  const size_t nckz= nc/2 + 1;
  const double dk= 2.0*M_PI/boxsize;
  const double knq= nc*M_PI/boxsize; // Nyquist frequency
  const double fac= pow(2*M_PI/boxsize, 1.5);
  const double fac_2pi3= 1.0/(8.0*M_PI*M_PI*M_PI);
  
  gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, seed);
  set_seedtable(nc, random_generator, seedtable);

  
  // clean the delta_k grid
  for(size_t ix=0; ix<local_nx; ix++)
   for(size_t iy=0; iy<nc; iy++)
    for(size_t iz=0; iz<nckz; iz++)
      for(int i=0; i<3; i++) {
	size_t index= (ix*nc + iy)*nckz + iz;
	psi_k[i][index][0] = 0;
	psi_k[i][index][1] = 0;
      }

  double kvec[3];
  for(size_t ix=0; ix<nc; ix++) {
    size_t iix = nc - ix;
    if(iix == nc)
      iix = 0;

    if(!((local_ix0 <= ix  && ix  < (local_ix0 + local_nx)) ||
	 (local_ix0 <= iix && iix < (local_ix0 + local_nx))))
      continue;
    
    for(size_t iy=0; iy<nc; iy++) {
      gsl_rng_set(random_generator, seedtable[ix*nc + iy]);
      
      for(size_t iz=0; iz<nc/2; iz++) {
	double phase= gsl_rng_uniform(random_generator)*2*M_PI;

	double mag= 1.0;
	if(!fix_amplitude) {
	  double ampl= 1.0;
	  do
	    ampl = gsl_rng_uniform(random_generator);
	  while(ampl == 0.0);
	  mag = -log(ampl);
	}

	if(ix == nc/2 || iy == nc/2 || iz == nc/2)
	  continue;
	if(ix == 0 && iy == 0 && iz == 0)
	  continue;
	
	if(ix < nc/2)
	  kvec[0]= dk*ix;
	else
	  kvec[0]= -dk*(nc - ix);
	
	if(iy < nc/2)
	  kvec[1]= dk*iy;
	else
	  kvec[1]= -dk*(nc - iy);
	
	if(iz < nc/2)
	  kvec[2]= dk*iz;
	else
	  kvec[2]= -dk*(nc - iz);
	
	double kmag2 = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];

	
#ifdef SPHEREMODE
	// select a sphere in k-space
	double kmag = sqrt(kmag2);
	if(kmag > knq)
	  continue;
#else
	if(fabs(kvec[0]) > knq)
	  continue;
	if(fabs(kvec[1]) > knq)
	  continue;
	if(fabs(kvec[2]) > knq)
	  continue;
#endif

	size_t index0= ((ix - local_ix0)*nc + iy)*nckz + iz;
	double delta2= mag*fac_2pi3*Pk[index0];

	double delta_k_mag= fac*sqrt(delta2);
	// delta_k_mag -- |delta_k| extrapolated to a=1
	// Displacement is extrapolated to a=1

	if(iz > 0) {
	  if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
	    size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
	    for(int i=0; i<3; i++) {
	      psi_k[i][index][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
	      psi_k[i][index][1]=  kvec[i]/kmag2*delta_k_mag*cos(phase);
	    }
	  }
	}
	else { // k=0 plane needs special treatment
	  if(ix == 0) {
	    if(iy >= nc/2)
	      continue;
	    else {
	      if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
		size_t iiy= nc - iy; // note: j!=0 surely holds at this point
		size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;				size_t iindex= ((ix - local_ix0)*nc + iiy)*nckz + iz;
		
		for(int i=0; i<3; i++) {
		  psi_k[i][index][0]=  -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][index][1]=   kvec[i]/kmag2*delta_k_mag*cos(phase);
		  
		  psi_k[i][iindex][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][iindex][1]= -kvec[i]/kmag2*delta_k_mag*cos(phase);
		}
	      }
	    }
	  }
	  else { // here comes i!=0 : conjugate can be on other processor!
	    if(ix >= nc/2)
	      continue;
	    else {
	      iix = nc - ix;
	      if(iix == nc)
		iix = 0;
	      int iiy = nc - iy;
	      if(iiy == nc)
		iiy = 0;
	      
	      if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
		size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
		for(int i=0; i<3; i++) {
		  psi_k[i][index][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][index][1]=  kvec[i]/kmag2*delta_k_mag*cos(phase);
		}
	      }
	      
	      if(local_ix0 <= iix && iix < (local_ix0 + local_nx)) {
		size_t index= ((iix - local_ix0)*nc + iiy)*nckz + iz;
		for(int i=0; i<3; i++) {
		  psi_k[i][index][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][index][1]= -kvec[i]/kmag2*delta_k_mag*cos(phase);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Inverse FFT to real space
  for(int i=0; i<3; ++i) {
    fft_psi[i]->mode= fft_mode_k;
    fft_psi[i]->execute_inverse();
  }
  
  gsl_rng_free(random_generator);  
}


void lpt_zeldovich_displacement_cic(vector<Particle>& v, const Float a)
{
  // Prerequisite:
  //     Psi(x) in fft_psi computed by lpt_generate_psi_x
  // Input:
  //     Particles in unperturbed location q
  // Output:
  //     Particle positions and velocities purturbed by Zeldovich approximation
  //     x = q + Psi(x)
  //     v = Dv*Psi(x)

  assert(comm_n_nodes() == 1); // Not MPIzed
  Float* psi[]=  {fft_psi[0]->fx, fft_psi[1]->fx, fft_psi[2]->fx};

  //const size_t nczr= 2*(nc/2 + 1);
  //const Float dx= boxsize/nc;
  const Float dx_inv= nc/boxsize;

  //double nmesh3_inv= 1.0/pow((double)nc, 3.0);
  uint64_t id= (uint64_t) local_ix0*nc*nc + 1;

  const Float D1= cosmology_D_growth(a);
  const Float Dv= cosmology_Dv_growth(a, D1);

  long double sum2= 0.0;

  const size_t n= v.size();
  
  for(Particle& p : v) {
    int ix0[3], ix1[3];
    double w0[3], w1[3];

    for(int k=0; k<3; ++k) {
      int ix= (int) (p.x[k]*dx_inv);
      w1[k]= p.x[k] - ix; // CIC weight for the right point
      w0[k]= 1 - w1[k];   //                    left point
      
      ix0[k]= ix % nc;       // left grid point
      ix1[k]= (ix + 1) % nc; // right grid point
    }
    
    // Interpolated Psi at particle position
    // This is not MPIzed, Psi can be in other MPI node
    
    for(int k=0; k<3; ++k) {
      Float dis=
	  w0[0]*w0[1]*w0[2]*psi[k][i_x(ix0[0], ix0[1], ix0[2])]
	+ w1[0]*w0[1]*w0[2]*psi[k][i_x(ix1[0], ix0[1], ix0[2])]
	+ w0[0]*w1[1]*w0[2]*psi[k][i_x(ix0[0], ix1[1], ix0[2])]
	+ w1[0]*w1[1]*w0[2]*psi[k][i_x(ix1[0], ix1[1], ix0[2])]
	+ w0[0]*w0[1]*w1[2]*psi[k][i_x(ix0[0], ix0[1], ix1[2])]
	+ w1[0]*w0[1]*w1[2]*psi[k][i_x(ix1[0], ix0[1], ix1[2])]
	+ w0[0]*w1[1]*w1[2]*psi[k][i_x(ix0[0], ix1[1], ix1[2])]
	+ w1[0]*w1[1]*w1[2]*psi[k][i_x(ix1[0], ix1[1], ix1[2])];

       p.x[k] += D1*dis;
       p.v[k]  = Dv*dis;

       sum2 += dis*dis;
     }
     p.id= id++;
  }


  sum2 /= (3.0*n);
  double f= cosmology_f_growth_rate(a);
  
  msg_printf(msg_verbose, "Zeldovich displacements CIC calculated.\n");
  msg_printf(msg_info, "scale_factor  a = %.15e\n", a);
  msg_printf(msg_info, "growth_factor D = %.15e\n", D1);
  msg_printf(msg_info, "goroth_rate   f = %.15e\n", f);
  msg_printf(msg_info, "displacement 1D rms %e\n", D1*sqrt(sum2));
  msg_printf(msg_info, "velocity 1D rms %e\n", Dv*sqrt(sum2));
}

void lpt_zeldovich_displacement_ngp(vector<Particle>& v, const Float a)
{
  // Prerequisite:
  //     Psi(x) in fft_psi computed by lpt_generate_psi_x
  // Input:
  //     Particles in unperturbed location q
  // Output:
  //     Particle positions and velocities purturbed by Zeldovich approximation
  //     x = q + Psi(x)
  //     v = Dv*Psi(x)

  assert(comm_n_nodes() == 1); // Not MPIzed
  Float* psi[]=  {fft_psi[0]->fx, fft_psi[1]->fx, fft_psi[2]->fx};

  const Float dx_inv= nc/boxsize;
  uint64_t id= (uint64_t) local_ix0*nc*nc + 1;

  const Float D1= cosmology_D_growth(a);
  const Float Dv= cosmology_Dv_growth(a, D1);

  long double sum2= 0.0;

  const size_t n= v.size();
  
  for(Particle& p : v) {
    int ix[3];

    for(int k=0; k<3; ++k)
      ix[k]= (int) (p.x[k]*dx_inv);

    // This is not MPIzed, Psi can be in other MPI node
    for(int k=0; k<3; ++k) {
      Float dis= psi[k][i_x(ix[0], ix[1], ix[2])];

      p.x[k] += D1*dis;
      p.v[k]  = Dv*dis;

      sum2 += dis*dis;
    }
    p.id= id++;
  }

  sum2 /= (3.0*n);
  double f= cosmology_f_growth_rate(a);
  
  msg_printf(msg_verbose, "Zeldovich displacements NGP calculated.\n");
  msg_printf(msg_info, "scale_factor  a = %.15e\n", a);
  msg_printf(msg_info, "growth_factor D = %.15e\n", D1);
  msg_printf(msg_info, "goroth_rate   f = %.15e\n", f);
  msg_printf(msg_info, "displacement 1D rms %e\n", D1*sqrt(sum2));
  msg_printf(msg_info, "velocity 1D rms %e\n", Dv*sqrt(sum2));
}

void lpt_set_offset(Float offset_)
{
  offset= offset_;
}

namespace {
  
void set_seedtable(const int nc, gsl_rng* random_generator,
		   unsigned int* const stable)
{
  // from N-GenIC
  assert(stable);
  

  for(int i=0; i<nc/2; i++) {
    for(int j=0; j<i; j++)
      stable[i*nc + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[j*nc + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[(nc - 1 - i)*nc + j] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[(nc - 1 - j)*nc + i] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[i*nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[j*nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[(nc - 1 - i)*nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[(nc - 1 - j)*nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);
  }
}





}
