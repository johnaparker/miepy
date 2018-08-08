#include "main.hpp"
//#include <cmath>
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
double associated_legendre(int n, int m, double z, bool derivative) {
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
