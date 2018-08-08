#ifndef GUARD_main_h
#define GUARD_main_h
#include <complex>

std::complex<double> spherical_hn(int n, double z, bool derivative=false);
std::complex<double> spherical_hn_2(int n, double z, bool derivative=false);
std::complex<double> test(int n, double z, bool derivative=false);

#endif
