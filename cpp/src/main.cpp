#include "main.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <omp.h>

using std::complex;

int add(int i, int j) {
    return i + j + 2;
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
