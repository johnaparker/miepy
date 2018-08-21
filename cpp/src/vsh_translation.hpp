#ifndef GUARD_vsh_translation_h
#define GUARD_vsh_translation_h

#include <complex>
#include <functional>
#include <array>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "vec.hpp"

enum class vsh_mode {
    outgoing,
    ingoing,
    incident,
    interior
};

std::function<std::complex<double>(int, double, bool)> get_zn(vsh_mode mode);

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
