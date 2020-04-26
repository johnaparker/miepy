#include "forces.hpp"
#include <complex>
#include <gsl/gsl_const_mksa.h>
#include "indices.hpp"
#include <cmath>

vec3 force(const Ref<const ComplexVector>& p_scat, const Ref<const ComplexVector>& p_inc,
        double k, double eps_b, double mu_b) {

    (void) mu_b;    //TODO: mu_b not being used

    static double eps_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;

    std::complex<double> Fxy = 0;
    std::complex<double> Fz = 0;
    double Axy = PI/pow(k,2)*eps_0*eps_b;
    double Az = -2*Axy;

    int rmax = p_scat.size()/2;
    int lmax = rmax_to_lmax(rmax);

    Eigen::Map<const ComplexVector> p(p_scat.data(), rmax);
    Eigen::Map<const ComplexVector> q(p_scat.data() + rmax, rmax);
    Eigen::Map<const ComplexVector> pi(p_inc.data(), rmax);
    Eigen::Map<const ComplexVector> qi(p_inc.data() + rmax, rmax);

    for (int n = 1; n < lmax + 1; n++) {
        for (int m = -n; m < n+1; m++) {
            int r = pow(n,2) + n + m - 1;
            // Fxy, term 1/3
            if (m != n) {
                double factor = Axy*sqrt((n+m+1)*(n-m))/double(n*(n+1));
                int r1 = r + 1;
                Fxy += factor*(2.0*p[r]*conj(q[r1])
                          - p[r]*conj(qi[r1])
                          - pi[r]*conj(q[r1])
                          + 2.0*q[r]*conj(p[r1])
                          - q[r]*conj(pi[r1])
                          - qi[r]*conj(p[r1]));
            }

            // Fz, term 1/2
            double factor = Az*m/double(n*(n+1));
            Fz += factor*(2.0*p[r]*conj(q[r])
                     - p[r]*conj(qi[r])
                     - pi[r]*conj(q[r]));

            if (n < lmax) {
                // Fxy, term 2/3
                factor = -Axy*sqrt((n+m+2)*(n+m+1)*n*(n+2)/double((2*n+3)*(2*n+1)))/(n+1);
                int r1 = pow(n+1,2) + (n+1) + m;

                Fxy += factor*(2.0*p[r]*conj(p[r1])
                          - p[r]*conj(pi[r1])
                          - pi[r]*conj(p[r1])
                          + 2.0*q[r]*conj(q[r1])
                          - q[r]*conj(qi[r1])
                          - qi[r]*conj(q[r1]));

                // Fxy, term 3/3
                factor = Axy*sqrt((n-m+1)*(n-m+2)*n*(n+2)/double((2*n+3)*(2*n+1)))/(n+1);
                r1 = pow(n+1,2) + (n+1) + m - 2;
                Fxy += factor*(2.0*p[r1]*conj(p[r])
                          - p[r1]*conj(pi[r])
                          - pi[r1]*conj(p[r])
                          + 2.0*q[r1]*conj(q[r])
                          - q[r1]*conj(qi[r])
                          - qi[r1]*conj(q[r]));

                // Fz, term 2/2
                factor = Az*sqrt((n-m+1)*(n+m+1)*n*(n+2)/double((2*n+3)*(2*n+1)))/(n+1);
                r1 = pow(n+1,2) + (n+1) + m - 1;
                Fz += factor*(2.0*p[r1]*conj(p[r])
                          - p[r1]*conj(pi[r])
                          - pi[r1]*conj(p[r])
                          + 2.0*q[r1]*conj(q[r])
                          - q[r1]*conj(qi[r])
                          - qi[r1]*conj(q[r]));
            }
        }
    }
    
    vec3 F;
    F << real(Fxy), imag(Fxy), real(Fz);
    return F;
}

vec3 torque(const Ref<const ComplexVector>& p_scat, const Ref<const ComplexVector>& p_inc,
        double k, double eps_b, double mu_b) {

    (void) mu_b;    //TODO: mu_b not being used
    static double eps_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;

    int rmax = p_scat.size()/2;
    int lmax = rmax_to_lmax(rmax);
    double A = -2*PI/pow(k,3)*eps_0*eps_b;
    vec3 T = vec3::Zero();

    Eigen::Map<const ComplexVector> p(p_scat.data(), rmax);
    Eigen::Map<const ComplexVector> q(p_scat.data() + rmax, rmax);
    Eigen::Map<const ComplexVector> pi(p_inc.data(), rmax);
    Eigen::Map<const ComplexVector> qi(p_inc.data() + rmax, rmax);

    for (int n = 1; n < lmax + 1; n++) {
        for (int m = -n; m < n+1; m++) {
            int r = pow(n,2) + n + m - 1;
            if (m != n) {
                double factor = -A*sqrt((n-m)*(n+m+1));
                int r1 = r + 1;
                
                // Tx
                T[0] += factor*real(p[r]*conj(p[r1])
                                  + q[r]*conj(q[r1])
                                  - 0.5*(p[r1]*conj(pi[r])
                                       + p[r]*conj(pi[r1])
                                       + q[r1]*conj(qi[r])
                                       + q[r]*conj(qi[r1])));

                // Ty
                T[1] += factor*imag(p[r]*conj(p[r1])
                                  + q[r]*conj(q[r1])
                                  + 0.5*(p[r1]*conj(pi[r])
                                       - p[r]*conj(pi[r1])
                                       + q[r1]*conj(qi[r])
                                       - q[r]*conj(qi[r1])));
            }

            // Tz
            double factor = A*m;
            T[2] += factor*(norm(p[r]) + norm(q[r])
                          - real(p[r]*conj(pi[r])
                               + q[r]*conj(qi[r])));
        }
    }

    return T;
}




