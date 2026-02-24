#include "special.hpp"
#include <cmath>
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#ifdef MIEPY_VALIDATE_GAUNT
#include <cassert>
#include <iostream>
#endif

using std::complex;
using namespace std::complex_literals;

double factorial(double n) {
    return std::tgamma(n+1);
}

// if n>z, the iterative jn's may not converge
complex<double> spherical_jn(int n, complex<double> z, bool derivative) {
    if (!derivative) {
        //return complex<double>(gsl_sf_bessel_jl(n, z), gsl_sf_bessel_yl(n, z));
        complex<double> sin_z = sin(z);
        complex<double> cos_z = cos(z);

        complex<double> jn_2p = sin_z/z;
        complex<double> jn = jn_2p;

        if (n > 0) {
            complex<double> jn_p = jn_2p/z - cos_z/z;
            jn = jn_p;

            for (int i = 2; i <= n; i++) {
                jn = double(2*i - 1)/z*jn_p - jn_2p;
                jn_2p = jn_p;
                jn_p = jn;
            }
        }

        return jn;
    }
    else {
        return spherical_jn(n-1, z) - double(n+1)/z*spherical_jn(n, z);
    }
}

double spherical_yn(int n, double z, bool derivative) {
    if (derivative)
        return spherical_yn(n-1, z) - (n+1)/z*spherical_yn(n, z);

    return gsl_sf_bessel_yl(n, z);
}

complex<double> spherical_hn(int n, double z, bool derivative) {
    if (!derivative) {
        //return complex<double>(gsl_sf_bessel_jl(n, z), gsl_sf_bessel_yl(n, z));
        double sin_z = sin(z);
        double cos_z = cos(z);

        double jn_2p = sin_z/z;
        double jn = jn_2p;

        double yn_2p = -cos_z/z;
        double yn = yn_2p;

        if (n > 0) {
            double jn_p = jn_2p/z - cos_z/z;
            double yn_p = yn_2p/z - jn_2p;
            jn = jn_p;
            yn = yn_p;

            for (int i = 2; i <= n; i++) {
                jn = (2*i - 1)/z*jn_p - jn_2p;
                jn_2p = jn_p;
                jn_p = jn;

                yn = (2*i - 1)/z*yn_p - yn_2p;
                yn_2p = yn_p;
                yn_p = yn;
            }
        }

        return std::complex<double>(jn, yn);
    }
    else {
        return spherical_hn(n-1, z) - (n+1)/z*spherical_hn(n, z);
    }
}

ComplexArray spherical_hn_recursion(int nmax, double z) {
    //return complex<double>(gsl_sf_bessel_jl(n, z), gsl_sf_bessel_yl(n, z));
    double sin_z = sin(z);
    double cos_z = cos(z);
    Array jn(nmax+1);
    Array yn(nmax+1);

    jn(0) = sin_z/z;
    yn(0) = -cos_z/z;

    jn(1) = jn(0)/z + yn(0);
    yn(1) = yn(0)/z - jn(0);

    for (int i = 2; i <= nmax; i++) {
        jn(i) = (2*i - 1)/z*jn(i-1) - jn(i-2);

        yn(i) = (2*i - 1)/z*yn(i-1) - yn(i-2);
    }

    ComplexArray ret = jn + 1.0i*yn;
    return ret;
}


complex<double> spherical_hn_2(int n, double z, bool derivative) {
    return std::conj(spherical_hn(n, z, derivative));
}

double riccati_1(int n, double z, bool derivative) {
    double jn = gsl_sf_bessel_jl(n, z);
    if (!derivative) {
        return z*jn;
    }
    else {
        double jn_d = gsl_sf_bessel_jl(n-1, z) - (n+1)/z*jn;
        return z*jn_d + jn;
    }
}

double riccati_2(int n, double z, bool derivative) {
    double yn = gsl_sf_bessel_yl(n, z);
    if (!derivative) {
        return -z*yn;
    }
    else {
        double yn_d = gsl_sf_bessel_yl(n-1, z) - (n+1)/z*yn;
        return -z*yn_d - yn;
    }
}

double riccati_3(int n, double z, bool derivative) {
    return riccati_2(n, z, derivative) - riccati_1(n, z, derivative);
}

// create and return Eigen array instead
// better performance for direct evaluation (especially when n=1) compared to scipy/sympy
// it may be worth included direct n=2 (still better than recursive)
// z is always equal to cos(theta) -- simplifies the n=(1,2) equations
double associated_legendre(int n, int m, double z, bool derivative) {
    if (!derivative) {
        if (n == 1) {
            if (m == 0)
                return z;
            if (m == 1)
                return sqrt(1 - pow(z,2));
            if (m == -1)
                return -sqrt(1 - pow(z,2))/2;
        }
        auto size = gsl_sf_legendre_array_n(n);
        auto index = gsl_sf_legendre_array_index(n, abs(m));
        double *result = new double[size];
        gsl_sf_legendre_array(GSL_SF_LEGENDRE_NONE, n, z, result);
        double leg = result[index];
        delete[] result;

        if (m < 0) {
            //double factor = std::tgamma(n + m + 1)/std::tgamma(n - m + 1);
            double factor = pow(-1, m)*factorial(n+m)/factorial(n-m);
            return factor*leg;
        }
        else {
            return leg;
        }
    }
    // TODO: implement
    else {
        return 0;
    }
}

Array associated_legendre_recursion(int nmax, double z) {
    auto size = gsl_sf_legendre_array_n(nmax);
    double *result = new double[size];
    gsl_sf_legendre_array(GSL_SF_LEGENDRE_NONE, nmax, z, result);

    Array ret(nmax*(nmax+2) + 1);

    for (int n = 0; n <= nmax; n++) {
        for (int m = -n; m <= n; m++) {
            auto index = gsl_sf_legendre_array_index(n, abs(m));
            double leg = result[index];

            if (m < 0) {
                double factor = pow(-1, m)*factorial(n+m)/factorial(n-m);
                leg *= factor;
            }

            int idx = n*(n+2) - n + m;
            ret(idx) = leg;
        }
    }

    delete[] result;
    return ret;
}

// combine below into pi_tau_func that computes both and returns as a tuple
// take the sympy reduced expressions for n=1,2 and use for direct computation
//          best case: recurse upward after n=2
double pi_func(int n, int m, double theta) {
    if (theta == 0) {
        if (m == 1)
            return n*(n+1)/2.0;
        else if (m == -1)
            return 1/2.0;
        else
            return 0;
    }
    else if (theta == PI) {
        if (m == 1)
            return pow(-1, n+1)*n*(n+1)/2.0;
        else if (m == -1)
            return pow(-1, n+1)/2.0;
        else
            return 0;
    }
    else {
        double z = cos(theta);
        return m/sin(theta)*associated_legendre(n, m, z);
    }

}

double tau_func(int n, int m, double theta) {
    if (theta == 0) {
        if (m == 1)
            return n*(n+1)/2.0;
        else if (m == -1)
            return -1/2.0;
        else
            return 0;
    }
    else if (theta == PI) {
        if (m == 1)
            return pow(-1, n)*n*(n+1)/2.0;
        else if (m == -1)
            return -pow(-1, n)/2.0;
        else
            return 0;
    }
    else {
        double z = cos(theta);
        auto size = gsl_sf_legendre_array_n(n);
        auto index = gsl_sf_legendre_array_index(n, abs(m));
        double *result = new double[size];
        double *result_p = new double[size];
        gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_NONE, n, z, result, result_p);
        double leg_p = -sin(theta)*result_p[index];
        delete[] result;
        delete[] result_p;

        if (m < 0) {
            //double factor = std::tgamma(n + m + 1)/std::tgamma(n - m + 1);
            double factor = pow(-1, m)*factorial(n+m)/factorial(n-m);
            return factor*leg_p;
        }
        else {
            return leg_p;
        }
    }
}

double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3) {
    return gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3);
}

double a_func(int m, int n, int u, int v, int p) {
    double numerator   = factorial(n+m)*factorial(v+u)*factorial(p-m-u);
    double denominator = factorial(n-m)*factorial(v-u)*factorial(p+m+u);
    double factor = pow(-1, m+u)*(2*p+1)*sqrt(numerator/denominator);

    double w1 = wigner_3j(n, v, p, 0, 0, 0);
    double w2 = wigner_3j(n, v, p, m, u, -m-u);

    return factor*w1*w2;
}

double b_func(int m, int n, int u, int v, int p) {
    double numerator   = factorial(n+m)*factorial(v+u)*factorial(p-m-u+1);
    double denominator = factorial(n-m)*factorial(v-u)*factorial(p+m+u+1);
    double factor = pow(-1, m+u)*(2*p+3)*sqrt(numerator/denominator);

    double w1 = wigner_3j(n, v, p, 0, 0, 0);
    double w2 = wigner_3j(n, v, p+1, m, u, -m-u);

    return factor*w1*w2;
}

wigner_3j_batch_result wigner_3j_batch(int j2, int j3, int m1, int m2, int m3) {
    wigner_3j_batch_result result;

    if (m1 + m2 + m3 != 0) {
        result.jmin = 0;
        result.jmax = -1;
        return result;
    }

    result.jmin = std::max(abs(j2 - j3), abs(m1));
    result.jmax = j2 + j3;
    int size = result.jmax - result.jmin + 1;

    if (size <= 0) {
        return result;
    }

    result.values.resize(size);

    if (size <= 2) {
        for (int i = 0; i < size; i++) {
            result.values[i] = wigner_3j(result.jmin + i, j2, j3, m1, m2, m3);
        }
        return result;
    }

    // Seed with 2 GSL calls
    result.values[0] = wigner_3j(result.jmin, j2, j3, m1, m2, m3);
    result.values[1] = wigner_3j(result.jmin + 1, j2, j3, m1, m2, m3);

    // Forward recursion via Schulten-Gordon three-term recurrence
    // j*s(j+1)*f(j+1) + (2j+1)*t(j)*f(j) + (j+1)*s(j)*f(j-1) = 0
    // where s(j) = sqrt[(j^2-(j2-j3)^2) * ((j2+j3+1)^2-j^2) * (j^2-m1^2)]
    //       t(j) = -m1*(j2(j2+1)-j3(j3+1)) + (m3-m2)*j*(j+1)
    // => f(j+1) = [-(j+1)*s(j)*f(j-1) - (2j+1)*t(j)*f(j)] / [j*s(j+1)]
    int d23 = j2 - j3;
    int s23p1 = j2 + j3 + 1;
    double j2j3_term = double(j2) * (j2 + 1) - double(j3) * (j3 + 1);
    double m_diff = m3 - m2;

    for (int i = 1; i < size - 1; i++) {
        int j = result.jmin + i;
        double j_d = j;
        double jp1_d = j + 1.0;

        double sj_sq = (j_d * j_d - double(d23) * d23) * (double(s23p1) * s23p1 - j_d * j_d)
                     * (j_d * j_d - double(m1) * m1);
        double sj = sqrt(std::max(0.0, sj_sq));

        double sjp1_sq = (jp1_d * jp1_d - double(d23) * d23) * (double(s23p1) * s23p1 - jp1_d * jp1_d)
                       * (jp1_d * jp1_d - double(m1) * m1);
        double sjp1 = sqrt(std::max(0.0, sjp1_sq));

        double A_coeff = jp1_d * sj;
        double B_coeff = (2.0 * j + 1) * (-double(m1) * j2j3_term + m_diff * j_d * jp1_d);
        double C_coeff = j_d * sjp1;

        if (std::abs(C_coeff) < 1e-30) {
            result.values[i + 1] = wigner_3j(j + 1, j2, j3, m1, m2, m3);
        } else {
            result.values[i + 1] = (-A_coeff * result.values[i - 1] - B_coeff * result.values[i]) / C_coeff;
        }
    }

    return result;
}

gaunt_result gaunt_batch(int m, int n, int u, int v) {
    gaunt_result result;

    int qmax_A = std::min({n, v, (n + v - abs(m + u)) / 2});
    int qmax_B = std::min({n, v, (n + v + 1 - abs(m + u)) / 2});

    result.a_vals.resize(qmax_A + 1);
    result.b_vals.resize(qmax_B + 1);
    result.b_vals[0] = 0.0;

    // Batch compute the two Wigner 3j sequences using column permutation symmetry:
    // W3j(n,v,p, 0,0,0) = W3j(p,n,v, 0,0,0) — batch over p as first argument
    auto w3j_zero = wigner_3j_batch(n, v, 0, 0, 0);

    // W3j(n,v,p, m,u,-(m+u)) = W3j(p,n,v, -(m+u),m,u) — batch over p as first argument
    auto w3j_mu = wigner_3j_batch(n, v, -(m + u), m, u);

    // Precompute shared factorial terms
    double fac_nm = factorial(n + m);
    double fac_n_m = factorial(n - m);
    double fac_vu = factorial(v + u);
    double fac_v_u = factorial(v - u);
    double sign = pow(-1, m + u);

    // Compute a_func values
    for (int q = 0; q <= qmax_A; q++) {
        int p = n + v - 2 * q;
        double fac_ratio = fac_nm * fac_vu * factorial(p - m - u) / (fac_n_m * fac_v_u * factorial(p + m + u));
        double factor = sign * (2 * p + 1) * sqrt(fac_ratio);

        double w1 = w3j_zero.values[p - w3j_zero.jmin];
        double w2 = w3j_mu.values[p - w3j_mu.jmin];

        result.a_vals[q] = factor * w1 * w2;

#ifdef MIEPY_VALIDATE_GAUNT
        double ref = a_func(m, n, u, v, p);
        double abs_err = std::abs(result.a_vals[q] - ref);
        double max_mag = std::max(std::abs(ref), std::abs(result.a_vals[q]));
        double rel_err = (max_mag > 1e-12) ? abs_err / max_mag : 0.0;
        if (rel_err > 1e-8) {
            std::cerr << "gaunt_batch a_vals mismatch: m=" << m << " n=" << n << " u=" << u << " v=" << v
                      << " p=" << p << " batch=" << result.a_vals[q] << " ref=" << ref
                      << " abs_err=" << abs_err << " rel_err=" << rel_err << "\n";
            assert(false);
        }
#endif
    }

    // Compute b_func values
    for (int q = 1; q <= qmax_B; q++) {
        int p = n + v - 2 * q;
        double fac_ratio = fac_nm * fac_vu * factorial(p - m - u + 1) / (fac_n_m * fac_v_u * factorial(p + m + u + 1));
        double factor = sign * (2 * p + 3) * sqrt(fac_ratio);

        double w1 = w3j_zero.values[p - w3j_zero.jmin];
        double w2 = w3j_mu.values[(p + 1) - w3j_mu.jmin];

        result.b_vals[q] = factor * w1 * w2;

#ifdef MIEPY_VALIDATE_GAUNT
        double ref = b_func(m, n, u, v, p);
        double abs_err = std::abs(result.b_vals[q] - ref);
        double max_mag = std::max(std::abs(ref), std::abs(result.b_vals[q]));
        double rel_err = (max_mag > 1e-12) ? abs_err / max_mag : 0.0;
        if (rel_err > 1e-8) {
            std::cerr << "gaunt_batch b_vals mismatch: m=" << m << " n=" << n << " u=" << u << " v=" << v
                      << " p=" << p << " batch=" << result.b_vals[q] << " ref=" << ref
                      << " abs_err=" << abs_err << " rel_err=" << rel_err << "\n";
            assert(false);
        }
#endif
    }

    return result;
}
