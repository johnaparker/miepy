#define NOMINMAX
#include "flux.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_particle_cross_sections(py::module &m) {
    // flux submodule
    m.def("particle_cross_sections", particle_cross_sections, 
           "p_scat"_a, "p_inc"_a, "p_src"_a, "k"_a, R"pbdoc(
        Returns (scattering, absorption, extinction) cross-sections of particle
    )pbdoc");
}

void bind_cluster_cross_sections(py::module &m) {
    m.def("cluster_cross_sections", cluster_cross_sections, 
           "p_cluster"_a, "p_src"_a, "k"_a, R"pbdoc(
        Returns (scattering, absorption, extinction) cross-sections of a cluster
    )pbdoc");
}
