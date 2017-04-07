#include <iostream>
#include <string>
#include <time.h>
#include "fs.h"


#include <gsl/gsl_rng.h>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

static void generate_random(const size_t np, const double boxsize,
			    unsigned int seed,
			    vector<Particle>& v);

static double* compute_pk_grid(const PowerSpectrum& ps,
			const size_t nc, const double boxsize);

int main(int argc, char* argv[])
{
  // Initialise MPI
  comm_mpi_init(&argc, &argv);

  if(comm_n_nodes() > 1) {
    cerr << "Error: only works with mpirun -n 1\n";
    return 1;
  }

  //
  // command-line options (Boost program_options)
  //
  options_description opt("options [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "power spectrum filename")
    ("nc", value<int>()->default_value(64), "number of grid per dimension")
    ("boxsize", value<double>()->default_value(1000.0), "length of box on a side")
    ("a", value<double>()->default_value(1.0, "1"), "scale factor")
    ("np", value<size_t>()->default_value(1000), "number of output particles")
    ("seed", value<unsigned int>()->default_value(0), "random seed")
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
  const double a= vm["a"].as<double>(); assert(a > 0);
  const size_t np= vm["np"].as<size_t>();
  const double boxsize= vm["boxsize"].as<double>(); assert(boxsize > 0.0);
  unsigned int seed= vm["seed"].as<unsigned int>();

  if(seed == 0)
    seed= (unsigned int) time(NULL);

  PowerSpectrum ps(filename.c_str());

  double* const pk_grid= compute_pk_grid(ps, nc, boxsize);

  lpt_init(nc, boxsize, 0);
  lpt_generate_psi_x(seed, pk_grid, false);

  vector<Particle> v;
  generate_random(np, boxsize, seed, v);

  lpt_zeldovich_displacement_ngp(v, a);


  // Finalize MPI
  comm_mpi_finalise();
}


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

  for(int ix=0; ix<nc; ++ix) {
   int kx= ix <= knq ? ix : ix - nc;
   double w_x= w(sin_fac*kx);
   for(int iy=0; iy<nc; ++iy) {
    int ky = iy <= knq ? iy : iy - nc;
    double w_xy= w_x*w(sin_fac*ky);
    for(int kz=0; kz<knq; ++kz) {
      double w_xyz= w_xy*w(sin_fac*kz);
      size_t index= (ix*nc + iy)*nckz + kz;
      double k= fac*sqrt(kx*kx + ky*ky + kz*kz);

      /*
      // Compute aliasing factor ???
      for(int ax=-na; ax<na; ++ax) {
       double kax= kx + ax*kp;
       for(int ay=-na; ay<na; ++ay) {
	double kay= ky + ay*kp;
	for(int az=-na; az<na; ++az) {
	  double kaz= kz + az*kp;
	  double ka= sqrt(kax*kax + kay*kay + kaz*kaz);
	  double kk= ka/k;

	      float w1= w(sin_fac*(kx+ax*kp))*
		        w(sin_fac*(ky+ay*kp))*
                        w(sin_fac*(kz+az*kp));
	      float w2= w1*w1;
	      float w4= w2*w2;
	      float Pkfac= pow(ka/k, neff);;
	      c2gg += w4*Pkfac;
	    }
	  }
	}
      */
      

      pkgrid[index] = ps.P(k)/w_xyz;
    }
   }
  }
  return pkgrid;
}
 
    
