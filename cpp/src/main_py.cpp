#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "main.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cpp, m) {
    m.doc() = R"pbdoc(
        C++ submodule of MiePy
        -----------------------

        .. currentmodule:: cpp

        .. autosummary::
           :toctree: _generate
    )pbdoc";

    py::module special = m.def_submodule("special", "special functions module");

    special.def("spherical_hn", py::vectorize(spherical_hn), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical hankel function of the first kind or its derivative
    )pbdoc");

    special.def("spherical_hn_2", py::vectorize(spherical_hn_2), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical hankel function of the second kind or its derivative
    )pbdoc");

    special.def("riccati_1", py::vectorize(riccati_1), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Riccati Bessel function of the first kind
    )pbdoc");

    special.def("riccati_2", py::vectorize(riccati_2), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Riccati Bessel function of the second kind
    )pbdoc");

    special.def("riccati_3", py::vectorize(riccati_3), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Riccati Bessel function of the third kind
    )pbdoc");

    special.def("associated_legendre", py::vectorize(associated_legendre), 
           "n"_a, "m"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Associated legendre function of integer order and degree
    )pbdoc");

    special.def("pi_func", py::vectorize(pi_func), 
           "n"_a, "m"_a, "theta"_a, R"pbdoc(
        pi special function that appears in the vector spherical harmonics
    )pbdoc");

    special.def("tau_func", py::vectorize(tau_func), 
           "n"_a, "m"_a, "theta"_a, R"pbdoc(
        tau special function that appears in the vector spherical harmonics
    )pbdoc");

    special.def("wigner_3j", wigner_3j, 
           "j1"_a, "j2"_a, "j3"_a, "m1"_a, "m2"_a, "m3"_a, R"pbdoc(
        Wigner 3-j coefficients
    )pbdoc");

    special.def("test", test, 
            py::arg(), py::arg(), py::arg("derivative") = false, R"pbdoc(
            Test function
    )pbdoc");

    special.def("test2", test2, 
            R"pbdoc(
            Test2 function
    )pbdoc");

}
