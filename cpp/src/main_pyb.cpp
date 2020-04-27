#define NOMINMAX
#include <pybind11/pybind11.h>

namespace py = pybind11;

// special submodule
void bind_spherical_jn(py::module &);
void bind_spherical_yn(py::module &);
void bind_spherical_hn(py::module &);
void bind_spherical_hn_2(py::module &);

void bind_riccati_1(py::module &);
void bind_riccati_2(py::module &);
void bind_riccati_3(py::module &);

void bind_associated_legendre(py::module &);
void bind_pi_func(py::module &);
void bind_tau_func(py::module &);

void bind_wigner_3j(py::module &);
void bind_a_func(py::module &);
void bind_b_func(py::module &);

// vsh functions submodule
void bind_enum_vsh_mode(py::module &);
void bind_Emn(py::module &);
void bind_vsh_electric(py::module &);
void bind_vsh_magnetic(py::module &);
void bind_expand_E_cluster(py::module &);

// vsh_translation submodule
void bind_vsh_translation(py::module &);
void bind_vsh_translation_numpy(py::module &);
void bind_vsh_translation_eigen(py::module &);
void bind_vsh_translation_lambda(py::module &);
void bind_vsh_translation_lambda_py(py::module &);

// interactions submodule
void bind_enum_solver(py::module &);
void bind_bicgstab(py::module &);
void bind_sphere_aggregate_tmatrix(py::module &);
void bind_particle_aggregate_tmatrix(py::module &);
void bind_reflection_matrix_nia(py::module &);
void bind_solve_linear_system(py::module &);

// forces submodule
void bind_force(py::module &);
void bind_torque(py::module &);

// flux submodule
void bind_particle_cross_sections(py::module &);
void bind_cluster_cross_sections(py::module &);

// decomposition submodule
void bind_trapz(py::module &);
void bind_trapz_2d(py::module &);
void bind_integrate_phase(py::module &);
void bind_grid_interpolate(py::module &);

// misc tests
void bind_test(py::module &);
void bind_test2(py::module &);
void bind_test3(py::module &);
void bind_combine_arrays(py::module &);

PYBIND11_MODULE(cpp, m) {
    m.doc() = R"pbdoc(
        C++ submodule of MiePy
        -----------------------

        .. currentmodule:: cpp

        .. autosummary::
           :toctree: _generate
    )pbdoc";

    // special submodule
    py::module special_m = m.def_submodule("special", "special functions module");

    bind_spherical_jn(special_m);
    bind_spherical_yn(special_m);
    bind_spherical_hn(special_m);
    bind_spherical_hn_2(special_m);

    bind_riccati_1(special_m);
    bind_riccati_2(special_m);
    bind_riccati_3(special_m);

    bind_associated_legendre(special_m);
    bind_pi_func(special_m);
    bind_tau_func(special_m);

    bind_wigner_3j(special_m);
    bind_a_func(special_m);
    bind_b_func(special_m);

    // vsh functions submodule
    py::module vsh_functions_m = m.def_submodule("vsh_functions", "vsh functions module");
    bind_enum_vsh_mode(vsh_functions_m);
    bind_Emn(vsh_functions_m);
    bind_vsh_electric(vsh_functions_m);
    bind_vsh_magnetic(vsh_functions_m);
    bind_expand_E_cluster(vsh_functions_m);

    // vsh_translation submodule
    py::module vsh_translation_m = m.def_submodule("vsh_translation", "vsh translation functions module");

    bind_vsh_translation(vsh_translation_m);
    bind_vsh_translation_numpy(vsh_translation_m);
    bind_vsh_translation_eigen(vsh_translation_m);
    bind_vsh_translation_lambda(vsh_translation_m);
    bind_vsh_translation_lambda_py(vsh_translation_m);

    // interactions submodule
    py::module interactions_m = m.def_submodule("interactions", "interactions functions module");

    bind_enum_solver(interactions_m);
    bind_bicgstab(interactions_m);
    bind_sphere_aggregate_tmatrix(interactions_m);
    bind_reflection_matrix_nia(interactions_m);
    bind_particle_aggregate_tmatrix(interactions_m);
    bind_solve_linear_system(interactions_m);

    // forces submodule
    py::module forces_m = m.def_submodule("forces", "force functions module");

    bind_force(forces_m);
    bind_torque(forces_m);

    // flux submodule
    py::module flux_m = m.def_submodule("flux", "flux functions module");

    bind_particle_cross_sections(flux_m);
    bind_cluster_cross_sections(flux_m);

    // decomposition submodule
    py::module decomposition_m = m.def_submodule("decomposition", "source decomposition module");

    bind_trapz(decomposition_m);
    bind_trapz_2d(decomposition_m);
    bind_integrate_phase(decomposition_m);
    bind_grid_interpolate(decomposition_m);

    // misc tests
    bind_test(special_m);
    bind_test2(special_m);
    bind_test3(vsh_translation_m);
    bind_combine_arrays(vsh_translation_m);
}
