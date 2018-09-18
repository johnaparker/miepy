#ifndef GUARD_vsh_translation_h
#define GUARD_vsh_translation_h

#include <complex>
#include <functional>
#include <array>
#include <map>
#include "vec.hpp"
#include "vsh_functions.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

struct vsh_cache {
    std::complex<double> factor;
    int qmax_A;
    int qmax_B;
    ComplexArray A;
    ComplexArray B;

    vsh_cache(int n, int m, int v, int u);
};

using vsh_cache_map = std::map<std::array<int,4>, vsh_cache>;

vsh_cache_map create_vsh_cache_map(int lmax);

void vsh_translation_insert_pair(Ref<ComplexMatrix> agg_tmatrix, const tmatrix_t& tmatrix, int i, int j, 
        double rad, double theta, double phi, double k, const vsh_cache_map& vsh_precompute);

void vsh_translation_insert_pair(Ref<ComplexMatrix> agg_tmatrix, const Ref<const ComplexMatrix>& mie, int i, int j, 
        double rad, double theta, double phi, double k, const vsh_cache_map& vsh_precompute);

py::array_t<double> combine_arrays(py::array_t<double> a, py::array_t<double> b);

std::array<std::complex<double>, 2> vsh_translation(
        int m, int n, int u, int v, double rad, double theta,
        double phi, double k, vsh_mode mode);

py::array_t<std::complex<double>> vsh_translation_numpy(
        int m, int n, int u, int v, py::array_t<double> rad, py::array_t<double> theta,
        py::array_t<double> phi, double k, vsh_mode mode);

AB_type vsh_translation_eigen(
        int m, int n, int u, int v, const Ref<const Array>& rad, const Ref<const Array>& theta,
        const Ref<const Array>& phi, double k, vsh_mode mode);

std::function<std::array<std::complex<double>, 2> (double, double, double, double)> vsh_translation_lambda(
        int m, int n, int u, int v, vsh_mode mode);

py::array_t<std::complex<double>> vsh_translation_lambda_py(
        int m, int n, int u, int v, py::array_t<double> rad, py::array_t<double> theta,
        py::array_t<double> phi, double k, vsh_mode mode);

double test3(int size, int cores);

#endif
