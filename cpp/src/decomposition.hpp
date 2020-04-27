#include "vec.hpp"
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <array>

std::complex<double> trapz(const Ref<const Array>& x, const Ref<const ComplexArray>& y);

std::complex<double> trapz_2d(const Ref<const Array>& x, const Ref<const Array>& y, const Ref<const ComplexMatrix>& f);

ComplexMatrix integrate_phase(const pybind11::array_t<double> rhat, const Ref<const vec3>& origin, double k, int rmax,
        const Ref<const Array>& theta, const Ref<const Array>& phi, const pybind11::array_t<std::complex<double>> p0);

pybind11::array_t<std::complex<double>> grid_interpolate(
        const std::array<pybind11::array_t<double>,2> grid,
        const pybind11::array_t<std::complex<double>> data, 
        const pybind11::array_t<double> pts,
        std::complex<double> fill_value = 0);
