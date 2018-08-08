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

    special.def("test", test, 
            py::arg(), py::arg(), py::arg("derivative") = false, R"pbdoc(
            Test function
    )pbdoc");

}
