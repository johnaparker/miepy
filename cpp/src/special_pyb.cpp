#define NOMINMAX
#include "special.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_spherical_jn(py::module &m) {
    m.def("spherical_jn", py::vectorize(spherical_jn), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical bessel function of the first kind or its derivative
    )pbdoc");
}

void bind_spherical_yn(py::module &m) {
    m.def("spherical_yn", py::vectorize(spherical_yn), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical bessel function of the second kind or its derivative
    )pbdoc");
}

void bind_spherical_hn(py::module &m) {
    m.def("spherical_hn", py::vectorize(spherical_hn), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical hankel function of the first kind or its derivative
    )pbdoc");
}

void bind_spherical_hn_2(py::module &m) {
    m.def("spherical_hn_2", py::vectorize(spherical_hn_2), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical hankel function of the second kind or its derivative
    )pbdoc");
}

void bind_riccati_1(py::module &m) {
    m.def("riccati_1", py::vectorize(riccati_1), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Riccati Bessel function of the first kind
    )pbdoc");
}

void bind_riccati_2(py::module &m) {
    m.def("riccati_2", py::vectorize(riccati_2), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Riccati Bessel function of the second kind
    )pbdoc");
}

void bind_riccati_3(py::module &m) {
    m.def("riccati_3", py::vectorize(riccati_3), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Riccati Bessel function of the third kind
    )pbdoc");
}

void bind_associated_legendre(py::module &m) {
    m.def("associated_legendre", py::vectorize(associated_legendre), 
           "n"_a, "m"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Associated legendre function of integer order and degree
    )pbdoc");
}

void bind_pi_func(py::module &m) {
    m.def("pi_func", py::vectorize(pi_func), 
           "n"_a, "m"_a, "theta"_a, R"pbdoc(
        pi special function that appears in the vector spherical harmonics
    )pbdoc");
}

void bind_tau_func(py::module &m) {
    m.def("tau_func", py::vectorize(tau_func), 
           "n"_a, "m"_a, "theta"_a, R"pbdoc(
        tau special function that appears in the vector spherical harmonics
    )pbdoc");
}

void bind_wigner_3j(py::module &m) {
    m.def("wigner_3j", wigner_3j, 
           "j1"_a, "j2"_a, "j3"_a, "m1"_a, "m2"_a, "m3"_a, R"pbdoc(
        Wigner 3-j coefficients
    )pbdoc");
}

void bind_a_func(py::module &m) {
    m.def("a_func", a_func, 
           "m"_a, "n"_a, "u"_a, "v"_a, "p"_a, R"pbdoc(
        a function (Gaunt coefficient)
    )pbdoc");
}

void bind_b_func(py::module &m) {
    m.def("b_func", b_func, 
           "m"_a, "n"_a, "u"_a, "v"_a, "p"_a, R"pbdoc(
        b function
    )pbdoc");
}

void bind_test(py::module &m) {
    m.def("test", test, 
            py::arg(), py::arg(), py::arg("derivative") = false, R"pbdoc(
            Test function
    )pbdoc");
}

void bind_test2(py::module &m) {
    m.def("test2", test2, 
            R"pbdoc(
            Test2 function
    )pbdoc");
}

