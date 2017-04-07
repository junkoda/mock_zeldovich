#ifndef CONST_H
#define CONST_H 1

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

// Internal unit is 1/h Mpc, km/s, solar mass.

namespace c {
constexpr double unit_length=   3.085678e24; // Mpc in cm
constexpr double unit_mass=     1.989e33;    // solar mass in g
constexpr double unit_velocity= 1.0e5;       // 1 km/s in cm/sec
constexpr double G_cgs=         6.672e-8;    // Grav const int cgs

constexpr double unit_time=     unit_length/unit_velocity;
constexpr double G=             G_cgs*unit_mass*(unit_time*unit_time)/
                                (unit_length*unit_length*unit_length);

constexpr double H0= 100.0;                  // km/s/(1/h Mpc)
constexpr double rho_crit_0= 3.0*H0*H0/(8.0*M_PI*G);

constexpr double delta_c= 1.686;

constexpr double cH0inv= 2997.92458;         // c/H0 [h^-1 Mpc]

}

#endif
