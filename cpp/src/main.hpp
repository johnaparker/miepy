#ifndef GUARD_main_h
#define GUARD_main_h
#include <complex>

//separate into 3 files: radial functions, angular functions, and a/b coefficients

std::complex<double> spherical_hn(int n, double z, bool derivative=false);
std::complex<double> spherical_hn_2(int n, double z, bool derivative=false);

double riccati_1(int n, double z, bool derivative=false);
double riccati_2(int n, double z, bool derivative=false);
double riccati_3(int n, double z, bool derivative=false);

double associated_legendre(int n, int m, double z, bool derivative=false);
double tau_func(int n, int m, double theta);
double pi_func(int n, int m, double theta);

double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3);

std::complex<double> test(int n, double z, bool derivative=false);
double test2();

#endif
