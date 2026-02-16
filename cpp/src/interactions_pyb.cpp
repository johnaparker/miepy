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

void bind_bicgstab_profiled(py::module &m) {
    m.def("bicgstab_profiled", [](const Ref<const ComplexMatrix>& A,
                const Ref<const ComplexVector>& b, int maxiter, double tolerance) {
                auto result = bicgstab_profiled(A, b, maxiter, tolerance);
                return py::make_tuple(result.solution, result.iterations, result.residual);
            },
        "A"_a, "b"_a, "maxiter"_a = 1000, "tolerance"_a = 1e-5, R"pbdoc(
        BiCGSTAB linear solver with profiling. Returns (solution, iterations, residual).
    )pbdoc");
}

void bind_bicgstab_preconditioned_profiled(py::module &m) {
    m.def("bicgstab_preconditioned_profiled", [](const Ref<const ComplexMatrix>& A,
                const Ref<const ComplexVector>& b,
                const Ref<const ComplexMatrix>& M_inv_blocks, int block_size,
                int maxiter, double tolerance) {
                auto result = bicgstab_preconditioned_profiled(A, b, M_inv_blocks, block_size, maxiter, tolerance);
                return py::make_tuple(result.solution, result.iterations, result.residual);
            },
        "A"_a, "b"_a, "M_inv_blocks"_a, "block_size"_a,
        "maxiter"_a = 1000, "tolerance"_a = 1e-5, R"pbdoc(
        Preconditioned BiCGSTAB solver with profiling. Returns (solution, iterations, residual).
    )pbdoc");
}

void bind_solve_linear_system(py::module &m) {
    m.def("solve_linear_system", [](const Ref<const ComplexMatrix>& agg_tmatrix,
                const Ref<const ComplexVector>& p_src, solver method,
                py::object M_inv_blocks_obj, int block_size) {
                if (M_inv_blocks_obj.is_none()) {
                    return solve_linear_system(agg_tmatrix, p_src, method);
                } else {
                    auto M_inv_blocks = M_inv_blocks_obj.cast<Ref<const ComplexMatrix>>();
                    return solve_linear_system_preconditioned(agg_tmatrix, p_src, M_inv_blocks, block_size, method);
                }
            },
        "agg_tmatrix"_a, "p_src"_a, "method"_a,
        "M_inv_blocks"_a = py::none(), "block_size"_a = 0, R"pbdoc(
        Solve the linear system:  p_inc = p_src - tmatrix*p_inc
        Optionally accepts block-diagonal preconditioner M_inv_blocks.
    )pbdoc");
}

void bind_reflection_matrix_nia(py::module &m) {
    m.def("reflection_matrix_nia", reflection_matrix_nia, 
           "positions"_a, "mie"_a, "k"_a, "reflection"_a, "z"_a, R"pbdoc(
        Obtain the reflection matrix for a cluster of spheres
    )pbdoc");
}
