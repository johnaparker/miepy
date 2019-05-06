#define NOMINMAX
#include "forces.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_force(py::module &m) {
    m.def("force", force, 
           "p_scat"_a, "p_inc"_a, "k"_a, 
           py::arg("eps_b") = 1, py::arg("mu_b") = 1, R"pbdoc(
        Force on particle from expansion coefficients
    )pbdoc");
}

void bind_torque(py::module &m) {
    m.def("torque", torque, 
           "p_scat"_a, "p_inc"_a, "k"_a, 
           py::arg("eps_b") = 1, py::arg("mu_b") = 1, R"pbdoc(
        Torque on particle from expansion coefficients
    )pbdoc");
}
