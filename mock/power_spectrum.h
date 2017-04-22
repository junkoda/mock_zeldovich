#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H 1

void compute_power_spectrum(const char filename[],
		    complex_t* fk,
		    const int nc, const double boxsize,
		    const double k_min, const double k_max, const double dk,
			    const int mas_correction);



#endif
