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
#include "power_spectrum.h"

using namespace std;

namespace {
  unsigned int* seedtable= 0;
  double boxsize;

  size_t nc= 0;
  size_t local_nx;
  size_t local_ix0;
  Float offset= 0.5;

  FFT* fft_psi[3];    // Zeldovichi displacement Psi_i
  FFT* fft_delta;

  void set_seedtable(const int nc, gsl_rng* random_generator,
		     unsigned int* const stable);
  void lpt_generate_psi_k(const unsigned long seed,
			  PowerSpectrum* const,
			  const int n_mas_correction,
			  const int fix_amplitude);

  void lpt_generate_delta_k(const unsigned long seed,
			    PowerSpectrum* const ps,
			    const int fix_amplitude);


  inline size_t idx(const int ix, const int iy, const int iz) {
    return 2*(nc/2 + 1)*(nc*ix + iy) + iz;
  }

}


void lpt_init(const int nc_, const double boxsize_, Mem* mem)
{
  // nc_: number of particles per dimension
  // boxsize_: box length on a side
  // Mem: Memory object for LPT, can be 0.
  //      If mem = 0, memory is allocated exclusively for LPT
  boxsize= boxsize_;
  nc= nc_;

  msg_printf(msg_debug, "lpt_init(nc= %d, boxsize= %.1lf)\n", nc, boxsize);

  if(mem)
    mem->use_from_zero(0);
  
  for(int i=0; i<3; i++)
    fft_psi[i]= new FFT("Psi_i", nc, mem, 0);

  fft_delta= new FFT("delta", nc, 0, 0);

  
  seedtable = (unsigned int *) malloc(nc*nc*sizeof(unsigned int));
  assert(seedtable);

  // checks
  local_nx= fft_psi[0]->local_nx;
  local_ix0= fft_psi[0]->local_ix0;
  
  for(int i=0; i<3; i++) {
    assert((size_t) fft_psi[i]->nc == nc);
    assert((size_t) fft_psi[i]->local_nx == local_nx);
    assert((size_t) fft_psi[i]->local_ix0 == local_ix0);
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

void lpt_test_power_spectrum(const char filename[],
			     const unsigned long seed,
			     PowerSpectrum* const ps,
			     const bool fix_amplitude)
{
  // Not MPI parallelised
  assert(comm_n_nodes() == 1);

  lpt_generate_delta_k(seed, ps, fix_amplitude);

  compute_power_spectrum(filename,
			 fft_psi[0]->fk, fft_psi[0]->nc,
			 boxsize, 0.0, 1.0, 0.01,
			 false);
}

void lpt_test_kpsi(const char filename[], 
		   const unsigned long seed,
		   PowerSpectrum* const ps,
		   const bool fix_amplitude,
		   const bool fft)
{
  // Test the accuracy of delta_k = ik.Psi
  
  // Not MPI parallelised
  assert(comm_n_nodes() == 1);
    

  if(nc == 0 || seedtable == 0) {
    msg_printf(msg_error,
	       "Error: lpt_init() not called before lpt_set_displacements");
    return;
  }
  
  msg_printf(msg_verbose, "TEST FFT\n");
  lpt_generate_psi_k(seed, ps, 0, fix_amplitude);

  // Convert Psi_k to realspace
  if(fft) {
    msg_printf(msg_verbose, "Fourier transforming Zeldovich displacements\n");
    for(int i=0; i<3; i++) {
      fft_psi[i]->execute_inverse();
    }

    msg_printf(msg_verbose, "FFT back to Fourier space\n");
    for(int i=0; i<3; i++) {
      fft_psi[i]->execute_forward();
    }
  }


  complex_t* psi[]=  {fft_psi[0]->fk, fft_psi[1]->fk, fft_psi[2]->fk};

  //const size_t nczr= 2*(nc/2 + 1);
  //const Float dx= boxsize/nc;
  const size_t ncz= nc/2 + 1;

  const double fac= 2.0*M_PI/boxsize;   // fundamental frequency
  double psi_fac= pow(boxsize, 1.5);
  if(fft)
    psi_fac /= (nc*nc*nc);

  // delta(k) = k.Psi (neglicting i in ik.Psi)
  for(int ix=0; ix<nc; ++ix) {
   double kx= ix <= nc/2 ? fac*ix : fac*(ix - (int) nc);

   for(int iy=0; iy<nc; ++iy) {
    double ky= iy <= nc/2 ? fac*iy : fac*(iy - (int) nc);
    
    //int iz0 = !(kx > 0.0 || (kx == 0.0 && ky > 0.0));

    for(int iz=0; iz<nc/2+1; ++iz) {
      double kz= fac*iz;
      size_t index= (ix*nc + iy)*ncz + iz;

      psi[0][index][0]  = psi_fac*kx*psi[0][index][0];
      psi[0][index][1]  = psi_fac*kx*psi[0][index][1];
      psi[0][index][0] += psi_fac*ky*psi[1][index][0];
      psi[0][index][1] += psi_fac*ky*psi[1][index][1];
      psi[0][index][0] += psi_fac*kz*psi[2][index][0];
      psi[0][index][1] += psi_fac*kz*psi[2][index][1];
    }
   }
  }

  compute_power_spectrum(filename,
			 fft_psi[0]->fk, nc,
			 boxsize, 0.0, 1.0, 0.01,
			 false);

}

void compute_zeldovich_linear(const char filename[], 
			      const unsigned long seed,
			      PowerSpectrum* const ps,
			      const double a,
			      const bool fix_amplitude)
{
  // Compute power spectrum of lienarised Zeldovich approximation
  
  // Not MPI parallelised
  assert(comm_n_nodes() == 1);
    

  if(nc == 0 || seedtable == 0) {
    msg_printf(msg_error,
	       "Error: lpt_init() not called before lpt_set_displacements");
    return;
  }
  
  msg_printf(msg_verbose, "Zeldovich linear\n");
  lpt_generate_psi_k(seed, ps, 0, fix_amplitude);

  // Convert Psi_k to realspace
  msg_printf(msg_verbose, "Fourier transforming Zeldovich displacements\n");
  for(int i=0; i<3; i++) {
    fft_psi[i]->execute_inverse();
  }

  const Float dx_inv= nc/boxsize;
  const size_t ncz= 2*(nc/2 + 1);

  Float* delta = fft_delta->fx;
  Float* psi[]=  {fft_psi[0]->fx, fft_psi[1]->fx, fft_psi[2]->fx};

  const Float D1= cosmology_D_growth(a);
  const double psi_fac= pow(boxsize, 1.5)/(nc*nc*nc);

  msg_printf(msg_info, "scale_factor  a = %.15le\n", a);
  msg_printf(msg_info, "growth factor D= %.15le\n", D1);


  size_t ix0[3], ix1[3];
  
  for(int ix=0; ix<nc; ++ix) {
    ix0[0]= ix;
    ix1[0]= (ix + 1) % nc;
    for(int iy=0; iy<nc; ++iy) {
      ix0[1]= iy;
      ix1[1]= (iy + 1) % nc;
      for(int iz=0; iz<nc; ++iz) {
	ix0[2]= iz;
	ix1[2]= (iz + 1) % nc;

	double d= 0.0;
	d +=   D1*psi[0][idx(ix0[0], ix0[1], ix0[2])]*dx_inv 
	     + D1*psi[1][idx(ix0[0], ix0[1], ix0[2])]*dx_inv 
	     + D1*psi[2][idx(ix0[0], ix0[1], ix0[2])]*dx_inv;

	d += - D1*psi[0][idx(ix1[0], ix0[1], ix0[2])]*dx_inv 
	     + D1*psi[1][idx(ix1[0], ix0[1], ix0[2])]*dx_inv 
	     + D1*psi[2][idx(ix1[0], ix0[1], ix0[2])]*dx_inv;

	d +=   D1*psi[0][idx(ix0[0], ix1[1], ix0[2])]*dx_inv
	     - D1*psi[1][idx(ix0[0], ix1[1], ix0[2])]*dx_inv
	     + D1*psi[2][idx(ix0[0], ix1[1], ix0[2])]*dx_inv;

	d += - D1*psi[0][idx(ix1[0], ix1[1], ix0[2])]*dx_inv
	     - D1*psi[1][idx(ix1[0], ix1[1], ix0[2])]*dx_inv
	     + D1*psi[2][idx(ix1[0], ix1[1], ix0[2])]*dx_inv;

	d +=   D1*psi[0][idx(ix0[0], ix0[1], ix1[2])]*dx_inv
	     + D1*psi[1][idx(ix0[0], ix0[1], ix1[2])]*dx_inv
	     - D1*psi[2][idx(ix0[0], ix0[1], ix1[2])]*dx_inv;

	d += - D1*psi[0][idx(ix1[0], ix0[1], ix1[2])]*dx_inv
	     + D1*psi[1][idx(ix1[0], ix0[1], ix1[2])]*dx_inv
	     - D1*psi[2][idx(ix1[0], ix0[1], ix1[2])]*dx_inv;

	d +=   D1*psi[0][idx(ix0[0], ix1[1], ix1[2])]*dx_inv
	     - D1*psi[1][idx(ix0[0], ix1[1], ix1[2])]*dx_inv
	     - D1*psi[2][idx(ix0[0], ix1[1], ix1[2])]*dx_inv;

	d += - D1*psi[0][idx(ix1[0], ix1[1], ix1[2])]*dx_inv
	     - D1*psi[1][idx(ix1[0], ix1[1], ix1[2])]*dx_inv
	     - D1*psi[2][idx(ix1[0], ix1[1], ix1[2])]*dx_inv;

	size_t grid_index= (ix *nc + iy )*ncz + iz;
	delta[grid_index] = psi_fac*d/4;

	/*
	size_t index=  (ix *nc + iy )*ncz + iz;

	size_t index0= (ix1*nc + iy )*ncz + iz;
	size_t index1= (ix *nc + iy1)*ncz + iz;
	size_t index2= (ix *nc + iy )*ncz + iz1;

	double d0= (psi[0][index] - psi[0][index0])*dx_inv;
	double d1= (psi[1][index] - psi[1][index1])*dx_inv;
	double d2= (psi[2][index] - psi[2][index2])*dx_inv;

	delta[index] = psi_fac*(d0 + d1 + d2);
	*/
      }
    }
  }

  fft_delta->mode= fft_mode_x;
	
  msg_printf(msg_verbose, "FFT delta to Fourier space\n");
  fft_delta->execute_forward();
  
  compute_power_spectrum(filename,
			 fft_delta->fk, nc,
			 boxsize, 0.0, 1.0, 0.01,
			 2);

  msg_printf(msg_verbose, "Linearlised Zeldovich power spectrum written\n");
}


void compute_zeldovich_nonlinear(const char filename[], 
			      const unsigned long seed,
			      PowerSpectrum* const ps,
				 const double a,
			      const bool fix_amplitude)
{
  // Compute power spectrum of lienarised Zeldovich approximation
  
  // Not MPI parallelised
  assert(comm_n_nodes() == 1);
    

  if(nc == 0 || seedtable == 0) {
    msg_printf(msg_error,
	       "Error: lpt_init() not called before lpt_set_displacements");
    return;
  }
  
  msg_printf(msg_verbose, "Zeldovich nonlinear\n");
  lpt_generate_psi_k(seed, ps, 0, fix_amplitude);

  // Convert Psi_k to realspace
  msg_printf(msg_verbose, "Fourier transforming Zeldovich displacements\n");
  for(int i=0; i<3; i++) {
    fft_psi[i]->execute_inverse();
  }

  const Float dx= boxsize/nc;
  const Float dx_inv= nc/boxsize;
  const size_t ncz= 2*(nc/2 + 1);
  size_t count= 0;

  Float* delta = fft_delta->fx;
  Float* psi[]=  {fft_psi[0]->fx, fft_psi[1]->fx, fft_psi[2]->fx};

  const Float D1= cosmology_D_growth(a);
  double psi_fac= pow(boxsize, 1.5)/(nc*nc*nc);

  msg_printf(msg_info, "scale_factor  a = %.15le\n", a);
  msg_printf(msg_info, "growth factor D= %.15le\n", D1);

  long double sum[]= {0.0, 0.0, 0.0};

  size_t ix0[3], ix1[3];
  
  for(int ix=0; ix<nc; ++ix) {
    ix0[0]= ix;
    ix1[0]= (ix + 1) % nc;
    for(int iy=0; iy<nc; ++iy) {
      ix0[1]= iy;
      ix1[1]= (iy + 1) % nc;
      for(int iz=0; iz<nc; ++iz) {
	ix0[2]= iz;
	ix1[2]= (iz + 1) % nc;

	/*
	double del[8][3];
	for(int ii=0; ii<2; ++ii) {
	  int ix0 = (ix + ii) % nc;
	  for(int jj=0; jj<2; ++jj) {
	    int iy0= (iy + jj) % nc;
	    for(int kk=0; kk<2; ++kk) {
	      int iz0= (iz + kk) % nc;
	      size_t psi_index= (ix0*nc + iy0)*nc + iz0;
	      del[ii + 2*jj + 4*kk][0]= D1*psi[0][psi_index]*dx_inv;
	      del[ii + 2*jj + 4*kk][1]= D1*psi[1][psi_index]*dx_inv;
	      del[ii + 2*jj + 4*kk][2]= D1*psi[2][psi_index]*dx_inv;
	    }
	  }
	}
	*/

	double d= 0.0;
	d += (0.5 + D1*psi[0][idx(ix0[0], ix0[1], ix0[2])]*dx_inv)*
	     (0.5 + D1*psi[1][idx(ix0[0], ix0[1], ix0[2])]*dx_inv)*
	     (0.5 + D1*psi[2][idx(ix0[0], ix0[1], ix0[2])]*dx_inv);

	d += (0.5 - D1*psi[0][idx(ix1[0], ix0[1], ix0[2])]*dx_inv)*
	     (0.5 + D1*psi[1][idx(ix1[0], ix0[1], ix0[2])]*dx_inv)*
	     (0.5 + D1*psi[2][idx(ix1[0], ix0[1], ix0[2])]*dx_inv);

	d += (0.5 + D1*psi[0][idx(ix0[0], ix1[1], ix0[2])]*dx_inv)*
	     (0.5 - D1*psi[1][idx(ix0[0], ix1[1], ix0[2])]*dx_inv)*
	     (0.5 + D1*psi[2][idx(ix0[0], ix1[1], ix0[2])]*dx_inv);

	d += (0.5 - D1*psi[0][idx(ix1[0], ix1[1], ix0[2])]*dx_inv)*
	     (0.5 - D1*psi[1][idx(ix1[0], ix1[1], ix0[2])]*dx_inv)*
	     (0.5 + D1*psi[2][idx(ix1[0], ix1[1], ix0[2])]*dx_inv);

	d += (0.5 + D1*psi[0][idx(ix0[0], ix0[1], ix1[2])]*dx_inv)*
	     (0.5 + D1*psi[1][idx(ix0[0], ix0[1], ix1[2])]*dx_inv)*
	     (0.5 - D1*psi[2][idx(ix0[0], ix0[1], ix1[2])]*dx_inv);

	d += (0.5 - D1*psi[0][idx(ix1[0], ix0[1], ix1[2])]*dx_inv)*
	     (0.5 + D1*psi[1][idx(ix1[0], ix0[1], ix1[2])]*dx_inv)*
	     (0.5 - D1*psi[2][idx(ix1[0], ix0[1], ix1[2])]*dx_inv);

	d += (0.5 + D1*psi[0][idx(ix0[0], ix1[1], ix1[2])]*dx_inv)*
	     (0.5 - D1*psi[1][idx(ix0[0], ix1[1], ix1[2])]*dx_inv)*
	     (0.5 - D1*psi[2][idx(ix0[0], ix1[1], ix1[2])]*dx_inv);

	d += (0.5 - D1*psi[0][idx(ix1[0], ix1[1], ix1[2])]*dx_inv)*
	     (0.5 - D1*psi[1][idx(ix1[0], ix1[1], ix1[2])]*dx_inv)*
	     (0.5 - D1*psi[2][idx(ix1[0], ix1[1], ix1[2])]*dx_inv);


	/*
	double d= 0.0;
	d += (0.5 + del[0][0])*(0.5 + del[0][1])*(0.5 + del[0][2]);
	d += (0.5 - del[1][0])*(0.5 + del[1][1])*(0.5 + del[1][2]);
	d += (0.5 + del[2][0])*(0.5 - del[2][1])*(0.5 + del[2][2]);
	d += (0.5 - del[3][0])*(0.5 - del[3][1])*(0.5 + del[3][2]);
	d += (0.5 + del[4][0])*(0.5 + del[4][1])*(0.5 - del[4][2]);
	d += (0.5 - del[5][0])*(0.5 + del[5][1])*(0.5 - del[5][2]);
	d += (0.5 + del[6][0])*(0.5 - del[6][1])*(0.5 - del[6][2]);
	d += (0.5 - del[7][0])*(0.5 - del[7][1])*(0.5 - del[7][2]);
	*/
	
	// grid index
	size_t grid_index= (ix *nc + iy )*ncz + iz;
	delta[grid_index] = psi_fac*(d - 1.0);

	// Statistics
	for(int k=0; k<3; ++k) {
	  double del= D1*psi[k][grid_index];
	  sum[k] += del*del;
	  if(abs(del) > dx)
	    count++;
	}
	
	//double d0= D1*(psi[0][index] - psi[0][index0])*dx_inv;
	//double d1= D1*(psi[1][index] - psi[1][index1])*dx_inv;
	//double d2= D1*(psi[2][index] - psi[2][index2])*dx_inv;

	//DEBUG!!!
	//delta[index]= psi_fac*(d0 + d1 + d2);
	//delta[index] = (1.0 + d0)*(1.0 + d1)*(1.0 + d2) - 1.0;
	//delta[index] = psi_fac*((1.0 + d0)*(1.0 + d1)*(1.0 + d2) - 1.0);

	fprintf(stderr, "%d %d %d %.15le %.15le %.15le\n",
		ix, iy, iz,
		psi[0][grid_index],
		psi[1][grid_index],
		psi[2][grid_index]);

      }
    }
  }

  for(int k=0; k<3; ++k) {
    double rms= sqrt(sum[k]/(nc*nc*nc));
    msg_printf(msg_info, "psi rms %e %e\n", rms, rms*dx_inv);
  }

  msg_printf(msg_info, "|psi| > dx %zu fraction %e\n",
	     count,
	     (double) count/(3*nc*nc*nc));



  fft_delta->mode= fft_mode_x;

  //char gfilename[64];
  //sprintf(gfilename, "grid_zeldovich_nonlinear_%zu.b", nc);
  //write_grid_real(gfilename, fft_delta);

  //sprintf(gfilename, "psi_zeldovich_nonlinear_%zu.b", nc);
  //write_psi_real(gfilename);


  
	
  msg_printf(msg_verbose, "FFT delta to Fourier space\n");
  fft_delta->execute_forward();
  
  compute_power_spectrum(filename,
			 fft_delta->fk, nc,
			 boxsize, 0.0, 1.0, 0.01,
			 2);

  msg_printf(msg_verbose, "Nonlinear Zeldovich power spectrum written\n");
}



void lpt_write_displacements(const char filename[], 
			   const unsigned long seed, PowerSpectrum* const ps,
			   const double a, const size_t np,
			   const bool fix_amplitude)
{
  // Place unperturbed particles in lattice
  // perturbe with Zeldovich approximation
  // subsample to np particles (expected value)
  // lattice initial condition
  if(nc == 0 || seedtable == 0) {
    msg_printf(msg_error,
	       "Error: lpt_init() not called before lpt_set_displacements");
    return;
  }
  
  msg_printf(msg_verbose, "Computing 2LPT\n");
  size_t np_local= local_nx*nc*nc;
  lpt_generate_psi_k(seed, ps, 0, fix_amplitude);

  // Convert Psi_k to realspace
  msg_printf(msg_verbose, "Fourier transforming 2LPT displacements\n");
  for(int i=0; i<3; i++) {
    fft_psi[i]->execute_inverse();
  }

  Float* psi[]=  {fft_psi[0]->fx, fft_psi[1]->fx, fft_psi[2]->fx};

  msg_printf(msg_verbose, "Setting particle grid and displacements\n");

  const size_t nczr= 2*(nc/2 + 1);
  const Float dx= boxsize/nc;

  //uint64_t id= (uint64_t) local_ix0*nc*nc + 1;

  //init rng
  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed + 10*comm_this_node());
  const double sample_rate = np > 0 ? np/pow((double) nc, 3.0) : 1.0;
  msg_printf(msg_debug, "sampling rate %e\n", sample_rate);

  // buffer
  vector<Particle> v;
  vector<Particle> receive;

  const int n_nodes= comm_n_nodes();
  cerr << "n_nodes= " << n_nodes << endl;
  vector<int> nrecv(n_nodes);
  vector<int> offsets(n_nodes);

  FILE* fp= 0;
  if(comm_this_node() == 0) {
    fp = fopen(filename, "w"); assert(fp);
  }
    
  
  const Float D1= cosmology_D_growth(a);
  //const Float Dv= cosmology_Dv_growth(a, D1);
  long double sum2 = 0.0;

  Particle p;
  Float x0[3];
  for(size_t ix=0; ix<local_nx; ix++) {
    v.clear();
    
    x0[0]= (local_ix0 + ix + offset)*dx;
    for(size_t iy=0; iy<nc; iy++) {
      x0[1]= (iy + offset)*dx;
      for(size_t iz=0; iz<nc; iz++) {
	x0[2]= (iz + offset)*dx;
	
	size_t index= (ix*nc + iy)*nczr + iz;
	for(int k=0; k<3; k++) {
	  Float dis=  psi[k][index];
	  
	  p.x[k]= x0[k] + D1*dis;
	  //p->v[k]= Dv*dis;	  
	}

	if(gsl_rng_uniform(rng) < sample_rate)
	  v.push_back(p);
      }
    }

    int nbuf= v.size();
    MPI_Gather(&nbuf, 1, MPI_INT,
	       nrecv.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    int sum= 0;
    if(comm_this_node() == 0) {
      for(int i=0; i<n_nodes; ++i) {
	offsets[i]= sizeof(Particle)*sum;
	sum += nrecv[i];
	//printf("receive from %d: %d %d\n", i, nrecv[i], offsets[i]);      
	nrecv[i] *= sizeof(Particle);
      }
      receive.resize(sum);
    }
    
    MPI_Gatherv(v.data(), nbuf*sizeof(Particle), MPI_BYTE,
		receive.data(), nrecv.data(), offsets.data(),
		MPI_BYTE, 0, MPI_COMM_WORLD);
    if(comm_this_node() == 0) {
      for(int i=0; i<sum; ++i)
	fprintf(fp, "%.15le %.15le %.15le\n",
		receive[i].x[0], receive[i].x[1], receive[i].x[2]);
    }
  }

  if(comm_this_node() == 0)
    fclose(fp);
  
  gsl_rng_free(rng);
    
  double f= cosmology_f_growth_rate(a);
  
  msg_printf(msg_info, "scale_factor  a = %.15le\n", a);
  msg_printf(msg_info, "growth factor D= %.15le\n", D1);
  msg_printf(msg_info, "growth_rate   f = %.15le\n", f);
  msg_printf(msg_info, "displacement_rms sigma = %.15lle\n",
	     sqrt(sum2/(3.0*np_local)));


  msg_printf(msg_verbose, "2LPT displacements calculated.\n");

  /*
  uint64_t nc64= nc;
  
  particles->np_local= np_local;
  particles->np_total= nc64*nc64*nc64;
  particles->a_x= a;
  particles->a_v= a;
  */

}

void lpt_set_displacements_ngp(const unsigned long seed,
			       PowerSpectrum* const ps,
			       const double a,
			       const size_t np,
			       const bool fix_amplitude,
			       Particles* particles)
{
  // random q and NGP interpolation for Psi(q)
  if(nc == 0 || seedtable == 0) {
    msg_printf(msg_error,
	       "Error: lpt_init() not called before lpt_set_displacements");
    return;
  }
  
  msg_printf(msg_verbose, "Computing 2LPT\n");
  assert(particles);
  size_t np_local= np;
  if(particles->np_allocated < np_local)
    msg_abort("Error: Not enough particles allocated to put initial particles\n"
	      "np_allocated= %lu < required %lu\n",
	      particles->np_allocated, np_local);
 
  lpt_generate_psi_k(seed, ps, 1, fix_amplitude);

  // Convert Psi_k to realspace
  msg_printf(msg_verbose, "Fourier transforming 2LPT displacements\n");
  for(int i=0; i<3; i++) {
    fft_psi[i]->execute_inverse();
  }

  Float* psi[]=  {fft_psi[0]->fx, fft_psi[1]->fx, fft_psi[2]->fx};

  msg_printf(msg_verbose, "Setting particle grid and displacements\n");

  const size_t ncz= 2*(nc/2 + 1);
  const Float dx_inv= nc/boxsize;
  Particle* p= particles->p;

  //uint64_t id= (uint64_t) local_ix0*nc*nc + 1;

  const Float D1= cosmology_D_growth(a);
  //const Float Dv= cosmology_Dv_growth(a, D1);
  long double sum2 = 0.0;
  
  assert(comm_n_nodes() == 1); // Not MPIzed
  
  //init rng
  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed + 1);

  for(size_t i=0; i<np; ++i) {
    Float x[3];
    int ix[3];

    for(int k=0; k<3; ++k) {
      x[k]= boxsize*gsl_rng_uniform(rng);
      ix[k]= (int) floor(x[k]*dx_inv);
      ix[k]= (ix[k] + nc) % nc;
    }
    
    size_t index= (ix[0]*nc + ix[1])*ncz + ix[2];
    for(int k=0; k<3; k++) {
      Float dis=  psi[k][index];
	  
      p->x[k]= x[k] + D1*dis;
      //p->v[k]= Dv*dis;
	  
      sum2 += dis*dis;	 
    }
    //p->id= id++;
	
    p++;
  }


  gsl_rng_free(rng);
    
  double f= cosmology_f_growth_rate(a);
  
  msg_printf(msg_info, "scale_factor  a = %.15le\n", a);
  msg_printf(msg_info, "growth factor D= %.15le\n", D1);
  msg_printf(msg_info, "growth_rate   f = %.15le\n", f);
  msg_printf(msg_info, "displacement_rms sigma = %.15lle\n",
	     sqrt(sum2/(3.0*np_local)));


  msg_printf(msg_verbose, "2LPT displacements calculated (random + NGP).\n");

  uint64_t nc64= nc;
  
  particles->np_local= np_local;
  particles->np_total= nc64*nc64*nc64;
  particles->a_x= a;
  particles->a_v= a;

}

void lpt_set_offset(Float offset_)
{
  offset= offset_;
  msg_printf(msg_info, "offset set to %e\n", offset);
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

/*
inline double w(const double tx)
{
  return tx == 0.0 ? 1.0 : sin(tx)/tx;
}
*/


void lpt_generate_psi_k(const unsigned long seed,
			PowerSpectrum* const ps,
			const int n_mas_correction,
			const int fix_amplitude)
{
  // Generates 1LPT (Zeldovich) displacements, Psi_k
  // from N-GenIC by Volker Springel
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

  //const double sin_fac= 0.5*boxsize/nc;

  
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

    if(ix < nc/2)
      kvec[0]= dk*ix;
    else
      kvec[0]= -dk*(nc - ix);

    if(!((local_ix0 <= ix  && ix  < (local_ix0 + local_nx)) ||
	 (local_ix0 <= iix && iix < (local_ix0 + local_nx))))
      continue;
    
    for(size_t iy=0; iy<nc; iy++) {
      gsl_rng_set(random_generator, seedtable[ix*nc + iy]);

      if(iy < nc/2)
	kvec[1]= dk*iy;
      else
	kvec[1]= -dk*(nc - iy);

      for(size_t iz=0; iz<nc/2; iz++) {
	double phase= gsl_rng_uniform(random_generator)*2*M_PI;
	double ampl;
	do
	  ampl = gsl_rng_uniform(random_generator);
	while(ampl == 0.0);


	if(ix == nc/2 || iy == nc/2 || iz == nc/2)
	  continue;
	if(ix == 0 && iy == 0 && iz == 0)
	  continue;
	
	if(iz < nc/2)
	  kvec[2]= dk*iz;
	else
	  kvec[2]= -dk*(nc - iz);

	double kmag2 = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];
	double kmag = sqrt(kmag2);
	
#ifdef SPHEREMODE
	// select a sphere in k-space
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

	double delta2;
	if(fix_amplitude)
	  delta2= fac_2pi3*ps->P(kmag);
	else
	  delta2= -log(ampl)*fac_2pi3*ps->P(kmag);
	
	double delta_k_mag= fac*sqrt(delta2); ///pow(w_xyz, n_mas_correction);
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
	      if((size_t) iiy == nc)
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

  for(int i=0; i<3; ++i) {
    fft_psi[i]->mode= fft_mode_k;
  }
  
  gsl_rng_free(random_generator);  
}


void lpt_generate_delta_k(const unsigned long seed,
			  PowerSpectrum* const ps,
			  const int fix_amplitude)
{
  // Generates delta(k)
  msg_printf(msg_verbose, "Generating delta_k...\n");
  msg_printf(msg_info, "Random Seed = %lu\n", seed);

  for(int i=0; i<3; i++) assert(fft_psi[i]);
  
  complex_t* psi_k[]= {fft_psi[0]->fk, fft_psi[1]->fk, fft_psi[2]->fk};

  const size_t nckz= nc/2 + 1;
  const double dk= 2.0*M_PI/boxsize;
  const double knq= nc*M_PI/boxsize; // Nyquist frequency
  //const double fac= pow(2*M_PI/boxsize, 1.5);
  //const double fac_2pi3= 1.0/(8.0*M_PI*M_PI*M_PI);
  
  gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, seed);
  set_seedtable(nc, random_generator, seedtable);

  
  // clean the delta_k grid
  for(size_t ix=0; ix<local_nx; ix++)
   for(size_t iy=0; iy<nc; iy++)
    for(size_t iz=0; iz<nckz; iz++) {
      size_t index= (ix*nc + iy)*nckz + iz;
      psi_k[0][index][0] = 0;
      psi_k[0][index][1] = 0;
    }

  double kvec[3];
  for(size_t ix=0; ix<nc; ix++) {
    size_t iix = nc - ix;
    if(iix == nc)
      iix = 0;

    if(ix < nc/2)
      kvec[0]= dk*ix;
    else
      kvec[0]= -dk*(nc - ix);

    if(!((local_ix0 <= ix  && ix  < (local_ix0 + local_nx)) ||
	 (local_ix0 <= iix && iix < (local_ix0 + local_nx))))
      continue;
    
    for(size_t iy=0; iy<nc; iy++) {
      gsl_rng_set(random_generator, seedtable[ix*nc + iy]);

      if(iy < nc/2)
	kvec[1]= dk*iy;
      else
	kvec[1]= -dk*(nc - iy);

      for(size_t iz=0; iz<nc/2; iz++) {
	double phase= gsl_rng_uniform(random_generator)*2*M_PI;
	double ampl;
	do
	  ampl = gsl_rng_uniform(random_generator);
	while(ampl == 0.0);


	if(ix == nc/2 || iy == nc/2 || iz == nc/2)
	  continue;
	if(ix == 0 && iy == 0 && iz == 0)
	  continue;
	
	if(iz < nc/2)
	  kvec[2]= dk*iz;
	else
	  kvec[2]= -dk*(nc - iz);

	double kmag2 = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];
	double kmag = sqrt(kmag2);
	
#ifdef SPHEREMODE
	// select a sphere in k-space
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

	double delta2;
	if(fix_amplitude)
	  delta2= ps->P(kmag);
	  //delta2= fac_2pi3*ps->P(kmag);
	else
	  delta2= -log(ampl)*ps->P(kmag);
	  //delta2= -log(ampl)*fac_2pi3*ps->P(kmag);
	
	//double delta_k_mag= fac*sqrt(delta2); ///pow(w_xyz, n_mas_correction);
	double delta_k_mag= sqrt(delta2); ///pow(w_xyz, n_mas_correction);
	// delta_k_mag -- |delta_k| extrapolated to a=1
	// Displacement is extrapolated to a=1

	if(iz > 0) {
	  if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
	    size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
	    psi_k[0][index][0]= delta_k_mag*cos(phase);
	    psi_k[0][index][1]= delta_k_mag*sin(phase);
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
		
		psi_k[0][index][0]=  delta_k_mag*cos(phase);
		psi_k[0][index][1]=  delta_k_mag*sin(phase);
		  
		psi_k[0][iindex][0]= delta_k_mag*cos(phase);
		psi_k[0][iindex][1]= delta_k_mag*sin(phase);
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
	      if((size_t) iiy == nc)
		iiy = 0;
	      
	      if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
		size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
		psi_k[0][index][0]= delta_k_mag*cos(phase);
		psi_k[0][index][1]= delta_k_mag*sin(phase);
	      }
	      
	      if(local_ix0 <= iix && iix < (local_ix0 + local_nx)) {
		size_t index= ((iix - local_ix0)*nc + iiy)*nckz + iz;
		psi_k[0][index][0]= delta_k_mag*cos(phase);
		psi_k[0][index][1]= delta_k_mag*sin(phase);
	      }
	    }
	  }
	}
      }
    }
  }

  fft_psi[0]->mode= fft_mode_k;
  
  gsl_rng_free(random_generator);  
}

/*
void lpt_compute_psi2_k(void)
{
  // Compute 2nd order Psi(2) from 1st order Psi
  //   Precondition Psi_k  in fft_psi[]->fk
  //   Result       Psi2_k in fft_psi2[]->fk (Fourier space)
  
  msg_printf(msg_verbose, "Computing 2LPT displacement fields...\n");

  const size_t nckz= nc/2 + 1;
  const double dk= 2.0*M_PI/boxsize;

  complex_t* psi_k[]= {fft_psi[0]->fk, fft_psi[1]->fk, fft_psi[2]->fk};

  //const double fac = pow(2*M_PI/boxsize, 1.5);
  double kvec[3];

  //
  // 2nd order LPT
  //
  complex_t* psi_ij_k[6];
  for(int i=0; i<6; i++) {
    //assert(fft_psi_ij[i]);
    psi_ij_k[i]= fft_psi_ij[i]->fk;
  }

  // Take derivative dPsi_i/dq_j in Fourier space
  for(size_t ix=0; ix<local_nx; ix++) {
    for(size_t iy=0; iy<nc; iy++) {
      for(size_t iz=0; iz<nckz; iz++) {
	size_t index= (ix*nc + iy)*nckz + iz;
	if((ix + local_ix0) < nc/2)
	  kvec[0]= dk*(ix + local_ix0);
	else
	  kvec[0]= -dk*(nc - (ix + local_ix0));
	      
	if(iy < nc/2)
	  kvec[1]= dk*iy;
	else
	  kvec[1]= -dk*(nc - iy);
	      
	if(iz < nc/2)
	  kvec[2]= dk*iz;
	else
	  kvec[2]= -dk*(nc - iz);
	      
	// Derivatives of ZA displacements
	// dPsi_i/dq_j -> sqrt(-1) k_j Psi_i(k)
	psi_ij_k[0][index][0]= -psi_k[0][index][1]*kvec[0]; // Psi_1,1
	psi_ij_k[0][index][1]=  psi_k[0][index][0]*kvec[0];

	psi_ij_k[1][index][0]= -psi_k[0][index][1]*kvec[1]; // Psi_1,2
	psi_ij_k[1][index][1]=  psi_k[0][index][0]*kvec[1];

	psi_ij_k[2][index][0]= -psi_k[0][index][1]*kvec[2]; // Psi_1,3
	psi_ij_k[2][index][1]=  psi_k[0][index][0]*kvec[2];
	      
	psi_ij_k[3][index][0]= -psi_k[1][index][1]*kvec[1]; // Psi_2,2
	psi_ij_k[3][index][1]=  psi_k[1][index][0]*kvec[1];

	psi_ij_k[4][index][0]= -psi_k[1][index][1]*kvec[2]; // Psi_2,3
	psi_ij_k[4][index][1]=  psi_k[1][index][0]*kvec[2];

	psi_ij_k[5][index][0]= -psi_k[2][index][1]*kvec[2]; // Psi_3,3
	psi_ij_k[5][index][1]=  psi_k[2][index][0]*kvec[2];
      }
    }
  }

  for(int i=0; i<6; i++) 
    fft_psi_ij[i]->mode= fft_mode_k;

  // Second-order displacement Psi(2)
  // div.Psi(2) = Sum_{i<j} [ Psi_i,j Psi_i,j - Psi_i,i Psi_j,j ]
  // in realspace

  msg_printf(msg_verbose, "Fourier transforming displacement gradient...\n");
  for(int i=0; i<6; i++) 
    fft_psi_ij[i]->execute_inverse();

  Float* psi_ij[6]; for(int i=0; i<6; i++) psi_ij[i]= fft_psi_ij[i]->fx;
  Float* const div_psi2= fft_div_psi2->fx; // == fft_psi_ij[3];

  size_t nczr= 2*(nc/2 + 1);
  for(size_t ix=0; ix<local_nx; ix++) {
    for(size_t iy=0; iy<nc; iy++) {
      for(size_t iz=0; iz<nc; iz++) {
	size_t index= (ix*nc + iy)*nczr + iz;

	div_psi2[index]=
	    psi_ij[0][index]*(psi_ij[3][index] + psi_ij[5][index])
	  + psi_ij[3][index]*psi_ij[5][index]
          - psi_ij[1][index]*psi_ij[1][index]
          - psi_ij[2][index]*psi_ij[2][index]
	  - psi_ij[4][index]*psi_ij[4][index];
      }
    }
  }

  fft_div_psi2->mode= fft_mode_x;

  // Solve Poisson eq. for div.Psi(2) in Fourier space
  msg_printf(msg_verbose, "Fourier transforming second order source...\n");
  
  fft_div_psi2->execute_forward();
  complex_t* div_psi2_k= fft_div_psi2->fk;
  complex_t* psi2_k[]= {fft_psi2[0]->fk, fft_psi2[1]->fk, fft_psi2[2]->fk};

  if(local_ix0 == 0) {
    for(int i=0; i<3; i++)
      psi2_k[i][0][0]= psi2_k[i][0][1] = 0.0;
    // avoid zero division kmag2 = 0
  }

  // Set 2nd-order Psi(2)
  for(size_t ix=0; ix<local_nx; ix++) {
    for(size_t iy=0; iy<nc; iy++) {
      int iz0= (ix + local_ix0 == 0) && (iy == 0); // skip kvec=(0,0,0)
      for(size_t iz=iz0; iz<nckz; iz++) {
	size_t index= (ix*nc + iy)*nckz + iz;
	if((ix + local_ix0) < nc/2)
	  kvec[0]=  dk*(ix + local_ix0);
	else
	  kvec[0]= -dk*(nc - (ix + local_ix0));
	
	if(iy < nc/2)
	  kvec[1]= dk*iy;
	else
	  kvec[1]= -dk*(nc - iy);
	
	if(iz < nc/2)
	  kvec[2] = dk*iz;
	else
	  kvec[2] = -dk*(nc - iz);
	
	double kmag2= kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];

#ifdef CHECK	
	assert(kmag2 > 0);
#endif	
	    
	// Psi(2)_k = div.Psi(2)_k * k / (sqrt(-1) k^2)
	for(int i=0; i<3; i++) {
	  psi2_k[i][index][0]=  div_psi2_k[index][1]*kvec[i]/kmag2;
	  psi2_k[i][index][1]= -div_psi2_k[index][0]*kvec[i]/kmag2;
	}
      }
    }
  }

  for(int i=0; i<3; i++) {
    fft_psi2[i]->mode= fft_mode_k;
  }

}
*/


}

void write_grid_real(const char filename[], FFT const * const fft)
{
  int nc= fft->nc;
  int ncz= 2*(nc/2 + 1);
  Float* p= fft->fx;
  
  FILE* fp= fopen(filename, "w");
  fwrite(&nc, sizeof(int), 1, fp);

  for(int ix=0; ix<nc; ix++) {
    for(int iy=0; iy<nc; ++iy) {
      fwrite(p, sizeof(Float), nc, fp);
      p += ncz;
    }
  }
  
  fwrite(&nc, sizeof(int), 1, fp);
  

  fclose(fp);
  msg_printf(msg_info, "Grid file %s written.\n", filename);
}

void write_psi_real(const char filename[])
{
  size_t ncz= 2*(nc/2 + 1);
  
  FILE* fp= fopen(filename, "w");
  fwrite(&nc, sizeof(int), 1, fp);

  for(int ix=0; ix<nc; ix++) {
    for(int iy=0; iy<nc; ++iy) {
      for(int iz=0; iz<nc; ++iz) {
	size_t index = (ix*nc + iy)*ncz;
	fwrite(&fft_psi[0]->fx[index], sizeof(Float), 1, fp);
	fwrite(&fft_psi[1]->fx[index], sizeof(Float), 1, fp);
	fwrite(&fft_psi[2]->fx[index], sizeof(Float), 1, fp);
      }
    }
  }
  
  fwrite(&nc, sizeof(int), 1, fp);
  

  fclose(fp);
  msg_printf(msg_info, "Psi file %s written.\n", filename);
}
