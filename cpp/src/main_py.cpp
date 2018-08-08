#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "main.hpp"
#include "other.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cpp, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: cpp

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("spherical_hn", py::vectorize(spherical_hn), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical hankel function of the first kind or its derivative
    )pbdoc");

    m.def("spherical_hn_2", py::vectorize(spherical_hn_2), 
           "n"_a, "z"_a, "derivative"_a=false, R"pbdoc(
        Spherical hankel function of the second kind or its derivative
    )pbdoc");

    m.def("test", test, 
            py::arg(), py::arg(), py::arg("derivative") = false, R"pbdoc(
            Test function
    )pbdoc");

    py::module m_sub = m.def_submodule("vsh", "A submodule of cpp");

    m_sub.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("multiply", &multiply, R"pbdoc(
        Multiply two numbers

        Some other explanation about the add function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
