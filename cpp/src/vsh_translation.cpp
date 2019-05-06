#define NOMINMAX
#include "vsh_translation.hpp"
#include "special.hpp"
#include "indices.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

using std::complex;
using namespace std::complex_literals;

vsh_cache::vsh_cache(int n, int m, int v, int u) {
    factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)));
    qmax_A = std::min({n, v, (n + v - abs(m+u))/2});
    qmax_B = std::min({n, v, (n + v + 1 - abs(m+u))/2});
    A = ComplexArray(qmax_A+1);
    B = ComplexArray(qmax_B+1);

    for (int q = 0; q < qmax_A+1; q++) {
        int p = n + v - 2*q;
        double aq = a_func(m, n, u, v, p);
        A(q) = aq*pow(1i, p)*double(n*(n+1) + v*(v+1) - p*(p+1));
    }

    B(0) = 0;
    for (int q = 1; q < qmax_B+1; q++) {
        int p = n + v - 2*q;
        double bq = b_func(m, n, u, v, p);
        B(q) = bq*pow(1i, p+1)*sqrt((pow(p+1,2) - pow(n-v,2))*(pow(n+v+1,2) - pow(p+1,2)));
    }
}

vsh_cache_map create_vsh_cache_map(int lmax) {
    vsh_cache_map ret;

    for (int n = 1; n < lmax+1; n++) {
        for (int m = -n; m < n+1; m++) {
            for (int v = 1; v < n+1; v++) {
                for (int u = -v; u < v+1; u++) {
                    std::array<int,4> key = {n, m, v, u};
                    ret.insert(vsh_cache_map::value_type(key, vsh_cache(n, m, v, u)));
                }
            }
        }
    }

    return ret;    
}

void vsh_translation_insert_pair(Ref<ComplexMatrix> agg_tmatrix, const tmatrix_t& tmatrix, int i, int j, 
        double rad, double theta, double phi, double k, const vsh_cache_map& vsh_precompute) {

    int rmax = tmatrix.dimensions()[1]/2;
    int lmax = rmax_to_lmax(rmax);

    int p_max = 2*lmax + 1;
    ComplexArray zn = spherical_hn_recursion(p_max, k*rad);
    Array Pnm = associated_legendre_recursion(p_max, cos(theta));

    ComplexMatrix A_matrix(agg_tmatrix.rows(), agg_tmatrix.cols());

    for (int n = 1; n < lmax+1; n++) {
        for (int m = -n; m < n+1; m++) {
            for (int v = 1; v < n+1; v++) {
                for (int u = -v; u < v+1; u++) {
                    if (n == v && u < -m)
                        continue;

                    m *= -1;

                    std::array<int,4> key = {n, m, v, u};
                    vsh_cache cache = vsh_precompute.at(key);
                    complex<double> factor = cache.factor;
                    int qmax_A = cache.qmax_A;
                    int qmax_B = cache.qmax_B;
                    ComplexArray A = cache.A;
                    ComplexArray B = cache.B;

                    complex<double> exp_phi = exp(1i*double(u+m)*phi);

                    complex<double> sum_term = 0;
                    for (int q = 0; q < qmax_A+1; q++) {
                        int p = n + v - 2*q;
                        int idx = p*(p+2) - p + (u+m);
                        sum_term += A(q)*Pnm(idx)*zn(p);
                    }

                    complex<double> A_translation = factor*exp_phi*sum_term;

                    sum_term = 0;
                    for (int q = 1; q < qmax_B+1; q++) {
                        int p = n + v - 2*q;
                        int idx = (p+1)*((p+1)+2) - (p+1) + (u+m);
                        sum_term += B(q)*Pnm(idx)*zn(p+1);
                    }

                    complex<double> B_translation = -factor*exp_phi*sum_term;
                    std::array<complex<double>,2> transfer{A_translation, B_translation};

                    m *= -1;

                    for (int a = 0; a < 2; a++) {
                        for (int b = 0; b < 2; b++) {
                            complex<double> val = transfer[(a+b)%2];
                            int idx = i*(2*rmax) + a*(rmax) + n*(n+2) - n + m - 1;
                            int idy = j*(2*rmax) + b*(rmax) + v*(v+2) - v + u - 1;
                            A_matrix(idx, idy) = val;

                            idx = j*(2*rmax) + a*(rmax) + n*(n+2) - n + m - 1;
                            idy = i*(2*rmax) + b*(rmax) + v*(v+2) - v + u - 1;
                            A_matrix(idx, idy) = pow(-1, n+v+a+b)*val;

                            if ((n == v && u != -m) || (n != v)) {
                                idx = i*(2*rmax) + b*(rmax) + v*(v+2) - v - u - 1;
                                idy = j*(2*rmax) + a*(rmax) + n*(n+2) - n - m - 1;
                                A_matrix(idx, idy) = pow(-1, m+u-a-b)*val;

                                idx = j*(2*rmax) + b*(rmax) + v*(v+2) - v - u - 1;
                                idy = i*(2*rmax) + a*(rmax) + n*(n+2) - n - m - 1;
                                A_matrix(idx, idy) = pow(-1, m+u+n+v)*val;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int n = 1; n < lmax+1; n++) {
    for (int m = -n; m < n+1; m++) {
    for (int a = 0; a < 2; a++) {
        for (int v = 1; v < lmax+1; v++) {
        for (int u = -v; u < v+1; u++) {
        for (int b = 0; b < 2; b++) {
            for (int vp = 1; vp < lmax+1; vp++) {
            for (int up = -vp; up < vp+1; up++) {
            for (int bp = 0; bp < 2; bp++) {

                int ida = a*(rmax) + n*(n+2) - n + m - 1;
                int idb = b*(rmax) + v*(v+2) - v + u - 1;
                int idc = bp*(rmax) + vp*(vp+2) - vp - up - 1;

                int idx = i*(2*rmax) + ida;
                int idy = j*(2*rmax) + idb;
                int mid = j*(2*rmax) + idc;

                agg_tmatrix(idx, idy) += A_matrix(idx, mid)*tmatrix(j, idc, idb);

                idx = j*(2*rmax) + ida;
                idy = i*(2*rmax) + idb;
                mid = i*(2*rmax) + idc;
                agg_tmatrix(idx, idy) += A_matrix(idx, mid)*tmatrix(i, idc, idb);
            }}}
        }}}
    }}}
}

void vsh_translation_insert_pair(Ref<ComplexMatrix> agg_tmatrix, const Ref<const ComplexMatrix>& mie, int i, int j, 
        double rad, double theta, double phi, double k, const vsh_cache_map& vsh_precompute) {

    int lmax = mie.cols()/2;
    int rmax = lmax_to_rmax(lmax);

    int p_max = 2*lmax + 1;
    ComplexArray zn = spherical_hn_recursion(p_max, k*rad);
    Array Pnm = associated_legendre_recursion(p_max, cos(theta));

    for (int n = 1; n < lmax+1; n++) {
        for (int m = -n; m < n+1; m++) {
            for (int v = 1; v < n+1; v++) {
                for (int u = -v; u < v+1; u++) {
                    if (n == v && u < -m)
                        continue;

                    m *= -1;

                    std::array<int,4> key = {n, m, v, u};
                    vsh_cache cache = vsh_precompute.at(key);
                    complex<double> factor = cache.factor;
                    int qmax_A = cache.qmax_A;
                    int qmax_B = cache.qmax_B;
                    ComplexArray A = cache.A;
                    ComplexArray B = cache.B;

                    complex<double> exp_phi = exp(1i*double(u+m)*phi);

                    complex<double> sum_term = 0;
                    for (int q = 0; q < qmax_A+1; q++) {
                        int p = n + v - 2*q;
                        int idx = p*(p+2) - p + (u+m);
                        sum_term += A(q)*Pnm(idx)*zn(p);
                    }

                    complex<double> A_translation = factor*exp_phi*sum_term;

                    sum_term = 0;
                    for (int q = 1; q < qmax_B+1; q++) {
                        int p = n + v - 2*q;
                        int idx = (p+1)*((p+1)+2) - (p+1) + (u+m);
                        sum_term += B(q)*Pnm(idx)*zn(p+1);
                    }

                    complex<double> B_translation = -factor*exp_phi*sum_term;
                    std::array<complex<double>,2> transfer{A_translation, B_translation};

                    m *= -1;

                    for (int a = 0; a < 2; a++) {
                        for (int b = 0; b < 2; b++) {
                            complex<double> val = transfer[(a+b)%2];
                            int idx = i*(2*rmax) + a*(rmax) + n*(n+2) - n + m - 1;
                            int idy = j*(2*rmax) + b*(rmax) + v*(v+2) - v + u - 1;
                            agg_tmatrix(idx, idy) = val*mie(j, b*lmax + v-1);

                            idx = j*(2*rmax) + a*(rmax) + n*(n+2) - n + m - 1;
                            idy = i*(2*rmax) + b*(rmax) + v*(v+2) - v + u - 1;
                            agg_tmatrix(idx, idy) = pow(-1, n+v+a+b)*val*mie(i, b*lmax + v-1);

                            if ((n == v && u != -m) || (n != v)) {
                                idx = i*(2*rmax) + b*(rmax) + v*(v+2) - v - u - 1;
                                idy = j*(2*rmax) + a*(rmax) + n*(n+2) - n - m - 1;
                                agg_tmatrix(idx, idy) = pow(-1, m+u-a-b)*val*mie(j, a*lmax + n-1);

                                idx = j*(2*rmax) + b*(rmax) + v*(v+2) - v - u - 1;
                                idy = i*(2*rmax) + a*(rmax) + n*(n+2) - n - m - 1;
                                agg_tmatrix(idx, idy) = pow(-1, m+u+n+v)*val*mie(i, a*lmax + n-1);
                            }
                        }
                    }
                }
            }
        }
    }
}


py::array_t<complex<double>> vsh_translation_lambda_py(
        int m, int n, int u, int v, py::array_t<double> rad, py::array_t<double> theta,
        py::array_t<double> phi, double k, vsh_mode mode) {


    auto fn = vsh_translation_lambda(m, n, u, v, mode);

    auto buf_rad = rad.unchecked<1>(), buf_theta = theta.unchecked<1>(), buf_phi = phi.unchecked<1>();
    int size = buf_theta.shape(0);

    py::array_t<complex<double>> result({2, size});
    auto buf = result.mutable_unchecked<2>();

    for (int i = 0; i < size; i++) {
        auto AB = fn(buf_rad(i), buf_theta(i), buf_phi(i), k);
        buf(0,i) = std::get<0>(AB);
        buf(1,i) = std::get<1>(AB);
    }

    return result;
}

std::function<std::array<complex<double>, 2> (double, double, double, double)> vsh_translation_lambda(
        int m, int n, int u, int v, vsh_mode mode) {
    
    m *= -1;
    auto zn = get_zn(mode);

    complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)));

    int qmax_A = std::min({n, v, (n + v - abs(m+u))/2});
    ComplexArray A(qmax_A+1);

    for (int q = 0; q < qmax_A+1; q++) {
        int p = n + v - 2*q;
        double aq = a_func(m, n, u, v, p);
        A(q) = aq*pow(1i, p)*double(n*(n+1) + v*(v+1) - p*(p+1));
    }


    int qmax_B = std::min({n, v, (n + v + 1 - abs(m+u))/2});
    ComplexArray B(qmax_B+1);
    B(0) = 0;

    for (int q = 1; q < qmax_B+1; q++) {
        int p = n + v - 2*q;
        double bq = b_func(m, n, u, v, p);
        B(q) = bq*pow(1i, p+1)*sqrt((pow(p+1,2) - pow(n-v,2))*(pow(n+v+1,2) - pow(p+1,2)));
    }


    return [m,n,u,v,zn,factor,qmax_A,qmax_B,A,B](double rad, double theta, double phi, double k) {
        double cos_theta = cos(theta);
        complex<double> exp_phi = exp(1i*double(u+m)*phi);
        complex<double> sum_term_A = 0, sum_term_B = 0;

        for (int q = 0; q < qmax_A+1; q++) {
            int p = n + v - 2*q;
            double Pnm_val = associated_legendre(p, u+m, cos_theta);
            sum_term_A += A(q)*Pnm_val*zn(p, k*rad, false);
        }
        
        for (int q = 1; q < qmax_B+1; q++) {
            int p = n + v - 2*q;
            double Pnm_val = associated_legendre(p+1, u+m, cos_theta);
            sum_term_B += B(q)*Pnm_val*zn(p+1, k*rad, false);
        }

        complex<double> A_translation = factor*exp_phi*sum_term_A;
        complex<double> B_translation = -factor*exp_phi*sum_term_B;

        return std::array<complex<double>,2>{A_translation, B_translation};
    };
}

// Instead, return tuple of ComplexArray
AB_type vsh_translation_eigen(
        int m, int n, int u, int v, const Ref<const Array>& rad, const Ref<const Array>& theta,
        const Ref<const Array>& phi, double k, vsh_mode mode) {

    size_t size = rad.size();
    AB_type result = AB_type::Zero(2, size);

    m *= -1;
    auto zn = get_zn(mode);

    double factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m));
    
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
            /(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)));
    

    py::array_t<complex<double>> exp_phi(buf_theta.shape(0));
    auto buf_exp_phi = exp_phi.mutable_unchecked<1>();

    py::array_t<double> cos_theta(buf_theta.shape(0));
    auto buf_cos_theta = cos_theta.mutable_unchecked<1>();

    for (ssize_t idx = 0; idx < buf_theta.shape(0); idx++) {
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

        for (ssize_t idx = 0; idx < buf_theta.shape(0); idx++) {
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

        for (ssize_t idx = 0; idx < buf_theta.shape(0); idx++) {
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

    for (ssize_t idx = 0; idx < buf_a.shape[0]; idx++) {
        buf(0,idx) = ptr_a[idx];
        buf(1,idx) = ptr_b[idx];
    }

    return result;
}

std::array<complex<double>, 2> vsh_translation(
        int m, int n, int u, int v, double rad, double theta,
        double phi, double k, vsh_mode mode) {
    
    m *= -1;
    auto zn = get_zn(mode);

    //complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            ///(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)))*exp(complex<double>(0, (u+m)*phi));
    complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)))*exp(1i*double(u+m)*phi);
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

    return std::array<complex<double>,2>{A_translation, B_translation};
}

double test3(int size, int cores) {
    //double sum = 0;

    //int n = 2, m = 2, u = 2, v = 2;
    //Array rad   = Array::Constant(1500, 0.3);
    //Array theta = Array::Constant(1500, 0.3);
    //Array phi   = Array::Constant(1500, 0.3);
    //double k = 1;
    //auto mode = vsh_mode::outgoing;

    //ComplexArray result = vsh_translation_eigen(m, n, u, v, rad, theta, phi, k, mode);
    //return sum;

    double sum = 0;

    AB_type result = AB_type::Zero(2, size);

    int n = 2, m = 2, u = 2, v = 2;
    double rad = 0.3, theta = 0.3, phi = 0.3;
    double k = 1;
    auto mode = vsh_mode::outgoing;

    auto fn = vsh_translation_lambda(m, n, u, v, mode);

    #pragma omp parallel for num_threads(cores)
    for (int i = 0; i <size; i++) {
        auto AB = fn(rad, theta, phi, k);
        result(0,i) = AB[0];
        result(1,i) = AB[1];
    }

    return sum;
}
