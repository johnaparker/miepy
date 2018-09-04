#ifndef GUARD_special_h
#define GUARD_special_h
#include <complex>
#include "vec.hpp"

//separate into 3 files: radial functions, angular functions, and a/b coefficients

double factorial(double n);

//double spherical_jn(int n, double z, bool derivative=false);
std::complex<double> spherical_jn(int n, std::complex<double> z, bool derivative=false);
double spherical_yn(int n, double z, bool derivative=false);
std::complex<double> spherical_hn(int n, double z, bool derivative=false);
std::complex<double> spherical_hn_2(int n, double z, bool derivative=false);
ComplexArray spherical_hn_recursion(int nmax, double z);

double riccati_1(int n, double z, bool derivative=false);
double riccati_2(int n, double z, bool derivative=false);
double riccati_3(int n, double z, bool derivative=false);

double associated_legendre(int n, int m, double z, bool derivative=false);
Array associated_legendre_recursion(int nmax, double z);
double tau_func(int n, int m, double theta);
double pi_func(int n, int m, double theta);

double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3);
double a_func(int m, int n, int u, int v, int p);
double b_func(int m, int n, int u, int v, int p);

std::complex<double> test(int n, double z, bool derivative=false);
double test2();

#endif
