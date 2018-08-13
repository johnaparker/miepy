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

ComplexMatrix vsh_translation_eigen(
        int m, int n, int u, int v, const Ref<const Array>& rad, const Ref<const Array>& theta,
        const Ref<const Array>& phi, double k, vsh_mode mode) {

    int size = rad.size();
    ComplexMatrix result = ComplexMatrix::Zero(2, size);

    m *= -1;
    auto zn = get_zn(mode);

    double factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /double(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)));
    
    ComplexArray exp_phi = exp(1i*double(u+m)*phi);
    Array cos_theta = cos(theta);

    int qmax = std::min({n, v, (n + v - abs(m+u))/2});

    for (int q = 0; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double aq = a_func(m, n, u, v, p);
        complex<double> A = aq*pow(1i, p)*double(n*(n+1) + v*(v+1) - p*(p+1));

        for (size_t idx = 0; idx < size; idx++) {
            double Pnm_val = associated_legendre(p, u+m, cos_theta(idx));
            complex<double> zn_val  = zn(p, k*rad(idx), false);
            result(0,idx) += factor*exp_phi(idx)*A*Pnm_val*zn_val;
        }
    }

    qmax = std::min({n, v, (n + v + 1 - abs(m+u))/2});

    for (int q = 1; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double bq = b_func(m, n, u, v, p);
        complex<double> A = bq*pow(1i, p+1)*sqrt((pow(p+1,2) - pow(n-v,2))*(pow(n+v+1,2) - pow(p+1,2)));

        for (size_t idx = 0; idx < size; idx++) {
            double Pnm_val = associated_legendre(p+1, u+m, cos_theta(idx));
            complex<double> zn_val  = zn(p+1, k*rad(idx), false);
            result(1,idx) -= factor*exp_phi(idx)*A*Pnm_val*zn_val;
        }
    }

    return result;
}

py::array_t<complex<double>> vsh_translation_numpy(
        int m, int n, int u, int v, py::array_t<double> rad, py::array_t<double> theta,
        py::array_t<double> phi, double k, vsh_mode mode) {

    auto buf_rad = rad.unchecked<1>(), buf_theta = theta.unchecked<1>(), buf_phi = phi.unchecked<1>();
    py::array_t<complex<double>> result({2, int(buf_theta.shape(0))});
    auto buf = result.mutable_unchecked<2>();


    m *= -1;
    auto zn = get_zn(mode);

    double factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /double(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)));
    

    py::array_t<complex<double>> exp_phi({int(buf_theta.shape(0))});
    auto buf_exp_phi = exp_phi.mutable_unchecked<1>();

    py::array_t<double> cos_theta({int(buf_theta.shape(0))});
    auto buf_cos_theta = cos_theta.mutable_unchecked<1>();

    for (size_t idx = 0; idx < buf_theta.shape(0); idx++) {
        buf_exp_phi(idx)   = exp(1i*double(u+m)*buf_phi(idx));
        buf_cos_theta(idx) = cos(buf_theta(idx));
        buf(0,idx) = 0;
        buf(1,idx) = 0;
    }

    int qmax = std::min({n, v, (n + v - abs(m+u))/2});

    for (int q = 0; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double aq = a_func(m, n, u, v, p);
        complex<double> A = aq*pow(1i, p)*double(n*(n+1) + v*(v+1) - p*(p+1));

        for (size_t idx = 0; idx < buf_theta.shape(0); idx++) {
            double Pnm_val = associated_legendre(p, u+m, buf_cos_theta(idx));
            complex<double> zn_val  = zn(p, k*buf_rad(idx), false);
            buf(0,idx) += factor*buf_exp_phi(idx)*A*Pnm_val*zn_val;
        }
    }

    qmax = std::min({n, v, (n + v + 1 - abs(m+u))/2});

    for (int q = 1; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double bq = b_func(m, n, u, v, p);
        complex<double> A = bq*pow(1i, p+1)*sqrt((pow(p+1,2) - pow(n-v,2))*(pow(n+v+1,2) - pow(p+1,2)));

        for (size_t idx = 0; idx < buf_theta.shape(0); idx++) {
            double Pnm_val = associated_legendre(p+1, u+m, buf_cos_theta(idx));
            complex<double> zn_val  = zn(p+1, k*buf_rad(idx), false);
            buf(1,idx) -= factor*buf_exp_phi(idx)*A*Pnm_val*zn_val;
        }
    }

    return result;
}

py::array_t<double> combine_arrays(py::array_t<double> a, py::array_t<double> b) {
    auto buf_a = a.request(), buf_b = b.request();

    if (buf_a.ndim != 1 || buf_b.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf_a.size != buf_b.size)
        throw std::runtime_error("Input shapes must match");

    py::array_t<double, py::array::c_style> result({2, int(buf_a.size)});
    auto buf = result.mutable_unchecked<2>();

    double *ptr_a = (double *) buf_a.ptr,
           *ptr_b = (double *) buf_b.ptr;

    for (size_t idx = 0; idx < buf_a.shape[0]; idx++) {
        buf(0,idx) = ptr_a[idx];
        buf(1,idx) = ptr_b[idx];
    }

    return result;
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

    int qmax = std::min({n, v, (n + v - abs(m+u))/2});
    complex<double> sum_term = 0;

    for (int q = 0; q < qmax+1; q++) {
        int p = n + v - 2*q;
        double aq = a_func(m, n, u, v, p);
        complex<double> A = aq*pow(1i, p)*double(n*(n+1) + v*(v+1) - p*(p+1));

        double Pnm_val = associated_legendre(p, u+m, cos_theta);
        sum_term += A*Pnm_val*zn(p, k*rad, false);
    }

    complex<double> A_translation = factor*sum_term;

    qmax = std::min({n, v, (n + v + 1 - abs(m+u))/2});
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
