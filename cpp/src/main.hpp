#ifndef GUARD_main_h
#define GUARD_main_h
#include <complex>

std::complex<double> spherical_hn(int n, double z, bool derivative=false);
std::complex<double> spherical_hn_2(int n, double z, bool derivative=false);
double associated_legendre(int n, int m, double z, bool derivative=false);
double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3);

std::complex<double> test(int n, double z, bool derivative=false);
double test2();

#endif
