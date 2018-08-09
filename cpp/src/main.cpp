#include "main.hpp"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <omp.h>

using std::complex;

int factorial(int n) {
    return (n == 0) ? 1 : factorial(n - 1) * n;
}

complex<double> spherical_hn(int n, double z, bool derivative) {
    if (!derivative) {
        return complex<double>(gsl_sf_bessel_jl(n, z), gsl_sf_bessel_yl(n, z));
    }
    else {
        return spherical_hn(n-1, z) - (n+1)/z*spherical_hn(n, z);
    }
}

complex<double> spherical_hn_2(int n, double z, bool derivative) {
    return std::conj(spherical_hn(n, z, derivative));
}

// create and return Eigen array instead
// better performance for direct evaluation (especially when n=1) compared to scipy/sympy
// it may be worth included direct n=2 (still better than recursive)
// z is always equal to cos(theta) -- simplifies the n=(1,2) equations
double associated_legendre(int n, int m, double z, bool derivative) {
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
        double factor = pow(-1, m)*factorial(n+m)/double(factorial(n-m));
        return factor*leg;
    }
    else {
        return leg;
    }
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
    else if (theta == M_PI) {
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
    else if (theta == M_PI) {
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

        if (m < 0) {
            //double factor = std::tgamma(n + m + 1)/std::tgamma(n - m + 1);
            double factor = pow(-1, m)*factorial(n+m)/double(factorial(n-m));
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

#pragma omp declare reduction( + : std::complex<double> : \
                       std::plus< std::complex<double> >( )) \
                       initializer(omp_priv = omp_orig)

complex<double> test(int n, double z, bool derivative) {
    complex<double> total = 0;
#pragma omp parallel for default(shared) reduction(+:total)
    for (int i = 0; i < 10000; i++) { 
        total += spherical_hn(n, z, derivative);
    }

    return total;
}

double test2() {
    double sum = 0;
    for (int i = 0; i < 10000; i++)
        sum += wigner_3j(2, 1, 3, 2, -1, -1);
    return sum;
}

//template<class Ret, class... Args>
//class Memorized 
//{
    //Map<tuple<Args...>, Ret> map;
    //function<Ret(Args...)> fn;
//public:
    //Memorized(function<Ret(Args...)> _fn) : fn(_fn) {}
    //Ret operator()(Args... args)
    //{
        //tuple = make_tuple(args...);
        //if(map contains tuple){
            //return map[tuple]
        //}
        //return map[tuple] = fn(args...);
    //}  
//} 
 
//auto fn = Memorized([](int x)->int { return x+1; });
