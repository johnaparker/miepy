#define NOMINMAX
#include "vsh_functions.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_enum_vsh_mode(py::module &m) {
    py::enum_<vsh_mode>(m, "vsh_mode")
        .value("outgoing", vsh_mode::outgoing)
        .value("ingoing",  vsh_mode::ingoing)
        .value("incident", vsh_mode::incident)
        .value("interior", vsh_mode::interior);
}

void bind_Emn(py::module &m) {
    m.def("Emn", Emn, "m"_a, "n"_a);
}

void bind_vsh_electric(py::module &m) {
    m.def("vsh_electric", vsh_electric, "n"_a, "m"_a,
            "mode"_a, "rad"_a, "theta"_a, "phi"_a, "k"_a, R"pbdoc(
        VSH electric (TM) mode
    )pbdoc");
}

void bind_vsh_magnetic(py::module &m) {
    m.def("vsh_magnetic", vsh_magnetic, "n"_a, "m"_a,
            "mode"_a, "rad"_a, "theta"_a, "phi"_a, "k"_a, R"pbdoc(
        VSH magnetic (TE) mode
    )pbdoc");
}

void bind_expand_E_cluster(py::module &m) {
    m.def("expand_E_cluster", expand_E_cluster,
            "pos"_a, "p_expand"_a, "mode"_a, "x"_a, "y"_a, "z"_a, "k"_a, R"pbdoc(
        Expand the electric field at a set of points from a cluster of particles
    )pbdoc");
}
