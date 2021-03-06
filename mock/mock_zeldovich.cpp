#include <iostream>
#include <string>
#include <time.h>
#include "fs.h"


#include <gsl/gsl_rng.h>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

/*
//static void generate_random(const size_t np, const double boxsize,
//			    unsigned int seed,
			    vector<Particle>& v);

static void generate_grid(const int nc, const double boxsize,
			  vector<Particle>& v);
			  
static double* compute_pk_grid(const PowerSpectrum& ps,
			const size_t nc, const double boxsize);
*/

int main(int argc, char* argv[])
{
  // Initialise MPI
  comm_mpi_init(&argc, &argv);

  //
  // command-line options (Boost program_options)
  //
  options_description opt("options [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "power spectrum filename")
    ("nc", value<int>()->default_value(64), "number of displacement grid per dimension")
    ("boxsize", value<double>()->default_value(1000.0), "length of box on a side")
    ("a", value<double>()->default_value(1.0, "1"), "scale factor")
    ("omega_m", value<double>()->default_value(0.308), "Omega_matter")
    ("np", value<double>()->default_value(0), "number of output particles")
    ("seed", value<unsigned int>()->default_value(0), "random seed. initialise by time if 0.")
    ("random", "place initial position random")
    ("output,o", value<string>(), "output ascii filename")
    ("fix-amplitude", "fix delta amplitude")
    ("lattice", "place unperturbed particles in lattice")
    ("random-ngp", "random unperturbed particles + NGP displacement interpolation")
    ("offset", value<double>()->default_value(0.5), "offset of initial lattice in units of grid spacing 0 <= offset < 1")
    ("test-delta", "compute power spectrum from delta(k)")
    ("test-kpsi", "compute power spectrum from delta(k) = ik.Psi(k)")
    ("test-fft", "compute power spectrum from delta(k) = ik.Psi(k) with FFT")
    ("zeldovich-linear", "compute power spectrum from Psi(x) with linear approximation")
    ("zeldovich-nonlinear", "compute power spectrum from Psi(x) with linear approximation")
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
  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  const double omega_m= vm["omega_m"].as<double>();
  const double a= vm["a"].as<double>(); assert(a > 0);
  const double boxsize= vm["boxsize"].as<double>(); assert(boxsize > 0.0);
  unsigned int seed= vm["seed"].as<unsigned int>();
  const bool fix_amplitude= vm.count("fix-amplitude");
  const size_t np= (size_t) vm["np"].as<double>();
    
  if(seed == 0) {
    seed= (unsigned int) time(NULL);
  }

  const double offset= vm["offset"].as<double>();
  lpt_set_offset(offset);

  //printf("seed = %u\n", seed);
  //printf("omega_m = %.4f\n", omega_m);


  PowerSpectrum ps(filename.c_str());

  //double* const pk_grid= compute_pk_grid(ps, nc, boxsize);

  cosmology_init(omega_m);
  lpt_init(nc, boxsize, 0);
 
  //Particles* particles;

  const string ofilename= vm["output"].as<string>();

  if(vm.count("lattice")) {
    lpt_write_displacements(ofilename.c_str(), seed, &ps, a, np, fix_amplitude);
  }
  else if(vm.count("random-ngp")) {
    abort();
    //lpt_set_displacements_ngp(seed, &ps, a, np, fix_amplitude);
  }
  else if(vm.count("test-delta")) {
    lpt_test_power_spectrum(ofilename.c_str(),
			    seed, &ps, fix_amplitude);
  }
  else if(vm.count("test-kpsi")) {
    lpt_test_kpsi(ofilename.c_str(),
		  seed, &ps, fix_amplitude, false);
  }
  else if(vm.count("test-fft")) {
    lpt_test_kpsi(ofilename.c_str(),
		  seed, &ps, fix_amplitude, true);
  }
  else if(vm.count("zeldovich-linear")) {
    compute_zeldovich_linear(ofilename.c_str(),
			     seed, &ps, a, fix_amplitude);
  }
  else if(vm.count("zeldovich-nonlinear")) {
    compute_zeldovich_nonlinear(ofilename.c_str(),
			     seed, &ps, a, fix_amplitude);
  }
  else {
    cerr << "Method not specified --lattice or --ngp";
    return 1;
  }
	


  // Finalize MPI
  comm_mpi_finalise();
}

/*
void generate_random(const size_t np, const double boxsize,
		     unsigned int seed,
		     vector<Particle>& v)
{
  v.reserve(np);
  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed + 1);

  Particle p;
  for(size_t i=0; i<np; ++i) {
    p.x[0]= boxsize*gsl_rng_uniform(rng);
    p.x[1]= boxsize*gsl_rng_uniform(rng);
    p.x[2]= boxsize*gsl_rng_uniform(rng);

    v.push_back(p);
  }
  gsl_rng_free(rng);
}

void generate_grid(const int nc, const double boxsize,
		   vector<Particle>& v)
{
  const double dx= boxsize/nc;
  const size_t np= (size_t) nc*nc*nc;
  v.reserve(np);

  Particle p;
  for(int ix=0; ix<nc; ++ix) {
    p.x[0]= ix*dx;
    for(int iy=0; iy<nc; ++iy) {
      p.x[1]= iy*dx;
      for(int iz=0; iz<nc; ++iz) {
	p.x[2]= iz*dx;
	v.push_back(p);
      }
    }
  }
}

// Interpolation window function
static inline double w(const double x)
{
  if(x == 0) return 1.0;

  double sinc = sin(x)/x;
  double sinc2 = sinc*sinc;
  
  return sinc2*sinc2;
}


double* compute_pk_grid(const PowerSpectrum& ps,
			const size_t nc, const double boxsize)
{
  // na: alias computed for (na <= n < na)^3 grids

  const int knq= nc/2;
  const size_t nckz= nc/2 + 1;
  double* pkgrid= (double*) malloc(sizeof(double)*nc*nc*nckz);
  double fac= 2.0*M_PI/boxsize;
  double sin_fac= M_PI/nc; 

  pkgrid[0] = 0.0;
  
  for(int ix=0; ix<nc; ++ix) {
   int kx= ix <= knq ? ix : ix - nc;
   double w_x= w(sin_fac*kx);
   for(int iy=0; iy<nc; ++iy) {
    int ky = iy <= knq ? iy : iy - nc;
    double w_xy= w_x*w(sin_fac*ky);

    int kz0 = (ix == 0 && iy == 0);
    for(int kz=kz0; kz<knq; ++kz) {
      double w_xyz= w_xy*w(sin_fac*kz);
      size_t index= (ix*nc + iy)*nckz + kz;
      double k= fac*sqrt(kx*kx + ky*ky + kz*kz);

      pkgrid[index] = ps.P(k); ///w_xyz;
    }
   }
  }
  return pkgrid;
}
*/ 
    
