#include "vsh_translation.hpp"
#include "main.hpp"
#include <cmath>
#include <algorithm>

using std::complex;
using namespace std::complex_literals;

std::function<complex<double>(int, double, bool)> get_zn(vsh_mode mode) {
    switch(mode) {
        case vsh_mode::outgoing:
            return spherical_hn;
        case vsh_mode::ingoing:
            return spherical_hn_2;
        case vsh_mode::incident:
        case vsh_mode::interior:
            return spherical_jn;

    }
}

std::tuple<complex<double>, complex<double>> vsh_translation(
        int m, int n, int u, int v, double rad, double theta,
        double phi, double k, vsh_mode mode) {
    
    m *= -1;
    auto zn = get_zn(mode);

    //complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            ///(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)))*exp(complex<double>(0, (u+m)*phi));
    complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /double(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)))*exp(1i*double(u+m)*phi);
    double cos_theta = cos(theta);

    int qmax = std::min({n, v, (n + v - abs(m+u)/2)});
    complex<double> sum_term = 0;

    for (int q = 0; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double aq = a_func(m, n, u, v, p);
        complex<double> A = aq*pow(1i, p)*double(n*(n+1) + v*(v+1) - p*(p+1));

        double Pnm_val = associated_legendre(p, u+m, cos_theta);
        sum_term += A*Pnm_val*zn(p, k*rad, false);
    }

    complex<double> A_translation = factor*sum_term;

    qmax = std::min({n, v, (n + v + 1 - abs(m+u)/2)});
    sum_term = 0;

    for (int q = 1; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double bq = b_func(m, n, u, v, p);
        complex<double> A = bq*pow(1i, p+1)*sqrt((pow(p+1,2) - pow(n-v,2))*(pow(n+v+1,2) - pow(p+1,2)));

        double Pnm_val = associated_legendre(p+1, u+m, cos_theta);
        sum_term += A*Pnm_val*zn(p+1, k*rad, false);
    }

    complex<double> B_translation = -factor*sum_term;

    return std::make_tuple(A_translation, B_translation);
}

double test3() {
    double sum = 0;

    int n = 2, m = 2, u = 2, v = 2;
    double rad = 0.3, theta = 0.3, phi = 0.3;
    double k =1;
    auto mode = vsh_mode::outgoing;

    m *= -1;
    auto zn = get_zn(mode);

    //complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            ///(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)))*exp(complex<double>(0, (u+m)*phi));
    complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /double(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)))*exp(1i*double(u+m)*phi);

    double cos_theta;
    for (int i = 0; i < 2; i++)
        cos_theta = cos(theta);

    int qmax = std::min({n, v, (n + v - abs(m+u)/2)});
    complex<double> sum_term = 0;

    for (int q = 0; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double aq = a_func(m, n, u, v, p);
        complex<double> A = aq*pow(1i, p)*double(n*(n+1) + v*(v+1) - p*(p+1));

        double Pnm_val;
        for (int i = 0; i < 2; i++)
            Pnm_val = associated_legendre(p, u+m, cos_theta);

        sum_term += A*Pnm_val*zn(p, k*rad, false);
    }

    complex<double> A_translation = factor*sum_term;

    qmax = std::min({n, v, (n + v + 1 - abs(m+u)/2)});
    sum_term = 0;

    for (int q = 1; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double bq = b_func(m, n, u, v, p);
        complex<double> A = bq*pow(1i, p+1)*sqrt((pow(p+1,2) - pow(n-v,2))*(pow(n+v+1,2) - pow(p+1,2)));

        double Pnm_val;
        for (int i = 0; i < 2; i++)
            Pnm_val = associated_legendre(p+1, u+m, cos_theta);

        sum_term += A*Pnm_val*zn(p+1, k*rad, false);
    }

    complex<double> B_translation = -factor*sum_term;

    return sum;
}
