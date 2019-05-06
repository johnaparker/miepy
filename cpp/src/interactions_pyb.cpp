#define NOMINMAX
#include "interactions.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_particle_aggregate_tmatrix(py::module &m) {
    m.def("particle_aggregate_tmatrix", [](const Ref<const position_t>& positions,
                Ref<ComplexMatrix> tmatrix, double k) {

                int Nparticles = tmatrix.rows();
                int cols = int(sqrt(tmatrix.cols()));
                const tmatrix_t tmatrix_map(tmatrix.data(), Nparticles, cols, cols);
                return particle_aggregate_tmatrix(positions, tmatrix_map, k);
            },
        "positions"_a, "tmatrix"_a, "k"_a, R"pbdoc(
        Obtain the particle-centered aggregate T-matrix for a cluster of particles
    )pbdoc");
}

void bind_sphere_aggregate_tmatrix(py::module &m) {
    m.def("sphere_aggregate_tmatrix", sphere_aggregate_tmatrix, 
           "positions"_a, "mie"_a, "k"_a, R"pbdoc(
        Obtain the particle-centered aggregate T-matrix for a cluster of spheres
    )pbdoc");
}

void bind_enum_solver(py::module &m) {
    py::enum_<solver>(m, "solver")
        .value("bicgstab", solver::bicgstab)
        .value("exact",  solver::exact);
}

void bind_bicgstab(py::module &m) {
    m.def("bicgstab", bicgstab, 
           "A"_a, "b"_a, "maxiter"_a, "tolerance"_a, R"pbdoc(
        BiCGSTAB linear solver
    )pbdoc");
}

void bind_solve_linear_system(py::module &m) {
    m.def("solve_linear_system", solve_linear_system, 
           "agg_tmatrix"_a, "p_src"_a, "method"_a, R"pbdoc(
        Solve the linear system:  p_inc = p_src - tmatrix*p_inc
    )pbdoc");
}

void bind_reflection_matrix_nia(py::module &m) {
    m.def("reflection_matrix_nia", reflection_matrix_nia, 
           "positions"_a, "mie"_a, "k"_a, "reflection"_a, "z"_a, R"pbdoc(
        Obtain the reflection matrix for a cluster of spheres
    )pbdoc");
}
