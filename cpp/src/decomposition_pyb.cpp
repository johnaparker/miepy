#include "decomposition.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_trapz(py::module &m) {
    m.def("trapz", trapz, 
           "x"_a, "y"_a, R"pbdoc(
        Numerically integrate 1D discrete samples using the trapezoidal rule
    )pbdoc");
}

void bind_trapz_2d(py::module &m) {
    m.def("trapz_2d", trapz_2d, 
           "x"_a, "y"_a, "f"_a, R"pbdoc(
        Numerically integrate 2D discrete samples using the trapezoidal rule
    )pbdoc");
}
