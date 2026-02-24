#define NOMINMAX
#include "vsh_translation.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;


void bind_vsh_translation(py::module &m) {
    m.def("vsh_translation", vsh_translation, 
           "m"_a, "n"_a, "u"_a, "v"_a, "rad"_a, "theta"_a,
           "phi"_a, "k"_a, "mode"_a, R"pbdoc(
        VSH translation coefficients
    )pbdoc");
}

void bind_vsh_translation_numpy(py::module &m) {
    m.def("vsh_translation_numpy", vsh_translation_numpy, 
           "m"_a, "n"_a, "u"_a, "v"_a, py::arg("rad").noconvert(),
           py::arg("theta").noconvert(), py::arg("phi").noconvert(),
           "k"_a, "mode"_a, R"pbdoc(
        VSH translation coefficients
    )pbdoc");
}

void bind_vsh_translation_eigen(py::module &m) {
    m.def("vsh_translation_eigen", vsh_translation_eigen, 
           "m"_a, "n"_a, "u"_a, "v"_a, "rad"_a, "theta"_a,
           "phi"_a, "k"_a, "mode"_a, R"pbdoc(
        VSH translation coefficients
    )pbdoc");
}

void bind_vsh_translation_lambda(py::module &m) {
    m.def("vsh_translation_lambda", vsh_translation_lambda, 
           "m"_a, "n"_a, "u"_a, "v"_a, "mode"_a, R"pbdoc(
        VSH translation coefficients
    )pbdoc");
}

void bind_vsh_translation_lambda_py(py::module &m) {
    m.def("vsh_translation_lambda_py", vsh_translation_lambda_py, 
           "m"_a, "n"_a, "u"_a, "v"_a, "rad"_a, "theta"_a,
           "phi"_a, "k"_a, "mode"_a, R"pbdoc(
        VSH translation coefficients
    )pbdoc");
}

void bind_create_vsh_cache_map(py::module &m) {
    m.def("create_vsh_cache_map", [](int lmax) {
        create_vsh_cache_map(lmax);
    }, "lmax"_a, R"pbdoc(
        Create VSH cache map for benchmarking purposes only.

        This is a timing-only wrapper: it calls the C++ create_vsh_cache_map()
        and discards the result. The returned cache map lives in C++ and cannot
        be passed back to Python. Use this to measure cache construction cost,
        not to obtain a usable cache object.
    )pbdoc");
}

