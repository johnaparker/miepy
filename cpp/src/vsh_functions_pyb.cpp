#include "vsh_functions.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_enum_vsh_mode(py::module &m) {
    py::enum_<vsh_mode>(m, "vsh_mode")
        .value("outgoing", vsh_mode::outgoing)
        .value("ingoing",  vsh_mode::ingoing)
        .value("incident", vsh_mode::incident)
        .value("interior", vsh_mode::interior);
}
