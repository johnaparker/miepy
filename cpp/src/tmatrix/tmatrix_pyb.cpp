#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "tmatrix/tmatrix_ebcm.hpp"

namespace py = pybind11;

void bind_tmatrix_submodule(py::module &m) {
    auto tmatrix_m = m.def_submodule("tmatrix", "EBCM T-matrix computation module");

    tmatrix_m.def("compute_spheroid",
        [](double axis_z, double axis_xy,
           double k, std::complex<double> n_rel,
           int lmax, int Nint,
           bool use_ds, bool complex_plane, double eps_z,
           bool extended_precision, bool conducting) -> py::array_t<std::complex<double>> {

            int rmax = lmax * (lmax + 2);
            std::vector<std::complex<double>> T_flat;

#if MIEPY_HAS_QUAD
            if (extended_precision) {
                tmatrix::SpheroidGeometry<__float128> geom(
                    static_cast<__float128>(axis_z), static_cast<__float128>(axis_xy));
                T_flat = tmatrix::compute_axisymmetric_tmatrix<__float128>(
                    geom, k, n_rel, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            } else
#endif
            {
                if (extended_precision) {
                    throw std::runtime_error(
                        "Extended precision (__float128) not available on this platform");
                }
                tmatrix::SpheroidGeometry<double> geom(axis_z, axis_xy);
                T_flat = tmatrix::compute_axisymmetric_tmatrix<double>(
                    geom, k, n_rel, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            }

            // Reshape to [2, rmax, 2, rmax]
            py::array_t<std::complex<double>> result({2, rmax, 2, rmax});
            auto buf = result.mutable_unchecked<4>();
            int idx = 0;
            for (int a1 = 0; a1 < 2; a1++)
                for (int r1 = 0; r1 < rmax; r1++)
                    for (int a2 = 0; a2 < 2; a2++)
                        for (int r2 = 0; r2 < rmax; r2++)
                            buf(a1, r1, a2, r2) = T_flat[idx++];

            return result;
        },
        py::arg("axis_z"), py::arg("axis_xy"),
        py::arg("k"), py::arg("n_rel"),
        py::arg("lmax"), py::arg("Nint") = 200,
        py::arg("use_ds") = true, py::arg("complex_plane") = false,
        py::arg("eps_z") = 0.95, py::arg("extended_precision") = false,
        py::arg("conducting") = false,
        R"pbdoc(
            Compute T-matrix for a spheroid using EBCM.

            Parameters
            ----------
            axis_z : float
                Semi-axis along symmetry axis (z).
            axis_xy : float
                Semi-axis perpendicular to symmetry axis.
            k : float
                Wavenumber in the medium.
            n_rel : complex
                Relative refractive index (sqrt(eps/eps_m)).
            lmax : int
                Maximum angular momentum order.
            Nint : int
                Number of quadrature points (default 200).
            use_ds : bool
                Use distributed sources (default True).
            complex_plane : bool
                Distribute sources in complex plane for oblate particles (default False).
            eps_z : float
                Source distribution extent parameter (default 0.95).
            extended_precision : bool
                Use __float128 precision internally (default False).
            conducting : bool
                If True, use PEC boundary conditions (default False).

            Returns
            -------
            T : ndarray, shape [2, rmax, 2, rmax]
                T-matrix in MiePy convention.
        )pbdoc"
    );

    tmatrix_m.def("compute_cylinder",
        [](double half_height, double radius,
           double k, std::complex<double> n_rel,
           int lmax, int Nint,
           bool use_ds, bool complex_plane, double eps_z,
           bool extended_precision, bool conducting) -> py::array_t<std::complex<double>> {

            int rmax = lmax * (lmax + 2);
            std::vector<std::complex<double>> T_flat;

#if MIEPY_HAS_QUAD
            if (extended_precision) {
                tmatrix::CylinderGeometry<__float128> geom(
                    static_cast<__float128>(half_height), static_cast<__float128>(radius));
                T_flat = tmatrix::compute_axisymmetric_tmatrix<__float128>(
                    geom, k, n_rel, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            } else
#endif
            {
                if (extended_precision) {
                    throw std::runtime_error(
                        "Extended precision (__float128) not available on this platform");
                }
                tmatrix::CylinderGeometry<double> geom(half_height, radius);
                T_flat = tmatrix::compute_axisymmetric_tmatrix<double>(
                    geom, k, n_rel, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            }

            py::array_t<std::complex<double>> result({2, rmax, 2, rmax});
            auto buf = result.mutable_unchecked<4>();
            int idx = 0;
            for (int a1 = 0; a1 < 2; a1++)
                for (int r1 = 0; r1 < rmax; r1++)
                    for (int a2 = 0; a2 < 2; a2++)
                        for (int r2 = 0; r2 < rmax; r2++)
                            buf(a1, r1, a2, r2) = T_flat[idx++];

            return result;
        },
        py::arg("half_height"), py::arg("radius"),
        py::arg("k"), py::arg("n_rel"),
        py::arg("lmax"), py::arg("Nint") = 200,
        py::arg("use_ds") = true, py::arg("complex_plane") = false,
        py::arg("eps_z") = 0.95, py::arg("extended_precision") = false,
        py::arg("conducting") = false,
        R"pbdoc(
            Compute T-matrix for a cylinder using EBCM.

            Parameters
            ----------
            half_height : float
                Half-height of the cylinder.
            radius : float
                Radius of the cylinder.
            k : float
                Wavenumber in the medium.
            n_rel : complex
                Relative refractive index.
            lmax : int
                Maximum angular momentum order.
            Nint : int
                Number of quadrature points (default 200).
            use_ds : bool
                Use distributed sources (default True).
            complex_plane : bool
                Distribute sources in complex plane (default False).
            eps_z : float
                Source distribution extent parameter (default 0.95).
            extended_precision : bool
                Use __float128 precision internally (default False).
            conducting : bool
                If True, use PEC boundary conditions (default False).

            Returns
            -------
            T : ndarray, shape [2, rmax, 2, rmax]
        )pbdoc"
    );

    // ---- Diagnostic Q-matrix functions ----
    tmatrix_m.def("diagnostic_Q_matrices_spheroid",
        [](double axis_z, double axis_xy,
           double k, std::complex<double> n_rel,
           int m, int lmax, int Nint,
           bool use_ds, bool complex_plane, double eps_z,
           bool extended_precision, bool conducting)
           -> std::pair<py::array_t<std::complex<double>>, py::array_t<std::complex<double>>> {

            int Nrank = lmax;
            int Nmax = (m == 0) ? Nrank : (Nrank - m + 1);
            if (Nmax <= 0) {
                throw std::runtime_error("m must be <= lmax; got m=" + std::to_string(m)
                                         + ", lmax=" + std::to_string(lmax));
            }
            int dim = 2 * Nmax;

            std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>> result;

#if MIEPY_HAS_QUAD
            if (extended_precision) {
                tmatrix::SpheroidGeometry<__float128> geom(
                    static_cast<__float128>(axis_z), static_cast<__float128>(axis_xy));
                result = tmatrix::diagnostic_Q_matrices_m<__float128>(
                    geom, k, n_rel, m, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            } else
#endif
            {
                if (extended_precision) {
                    throw std::runtime_error(
                        "Extended precision (__float128) not available on this platform");
                }
                tmatrix::SpheroidGeometry<double> geom(axis_z, axis_xy);
                result = tmatrix::diagnostic_Q_matrices_m<double>(
                    geom, k, n_rel, m, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            }

            // Determine Q31 and Q11 dimensions
            // DS mode: Q31 is [2*Nrank, 2*Nrank], Q11 is [2*Nmax, 2*Nrank]
            // Non-DS: both are [2*Nmax, 2*Nmax]
            int rows31, cols31, rows11, cols11;
            if (use_ds) {
                rows31 = 2 * Nrank; cols31 = 2 * Nrank;
                rows11 = 2 * Nmax;  cols11 = 2 * Nrank;
            } else {
                rows31 = 2 * Nmax; cols31 = 2 * Nmax;
                rows11 = 2 * Nmax; cols11 = 2 * Nmax;
            }

            py::array_t<std::complex<double>> Q31({rows31, cols31});
            py::array_t<std::complex<double>> Q11({rows11, cols11});
            auto buf31 = Q31.mutable_unchecked<2>();
            auto buf11 = Q11.mutable_unchecked<2>();
            int idx31 = 0;
            for (int i = 0; i < rows31; i++)
                for (int j = 0; j < cols31; j++)
                    buf31(i, j) = result.first[idx31++];
            int idx11 = 0;
            for (int i = 0; i < rows11; i++)
                for (int j = 0; j < cols11; j++)
                    buf11(i, j) = result.second[idx11++];

            return {Q31, Q11};
        },
        py::arg("axis_z"), py::arg("axis_xy"),
        py::arg("k"), py::arg("n_rel"),
        py::arg("m"), py::arg("lmax"), py::arg("Nint") = 200,
        py::arg("use_ds") = true, py::arg("complex_plane") = false,
        py::arg("eps_z") = 0.95, py::arg("extended_precision") = false,
        py::arg("conducting") = false,
        R"pbdoc(
            Diagnostic: return Q31 and Q11 matrices for a spheroid at azimuthal mode m.

            Uses the same code path as compute_spheroid but returns the intermediate
            Q matrices instead of the final T-matrix.

            Returns
            -------
            Q31, Q11 : tuple of ndarray, each shape [2*Nmax, 2*Nmax]
                Raw Q matrices for the given azimuthal mode.
        )pbdoc"
    );

    tmatrix_m.def("compute_rounded_cylinder",
        [](double half_height, double radius,
           double k, std::complex<double> n_rel,
           int lmax, int Nint,
           bool use_ds, bool complex_plane, double eps_z,
           bool extended_precision, bool conducting) -> py::array_t<std::complex<double>> {

            int rmax = lmax * (lmax + 2);
            std::vector<std::complex<double>> T_flat;

#if MIEPY_HAS_QUAD
            if (extended_precision) {
                tmatrix::RoundedCylinderGeometry<__float128> geom(
                    static_cast<__float128>(half_height), static_cast<__float128>(radius));
                T_flat = tmatrix::compute_axisymmetric_tmatrix<__float128>(
                    geom, k, n_rel, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            } else
#endif
            {
                if (extended_precision) {
                    throw std::runtime_error(
                        "Extended precision (__float128) not available on this platform");
                }
                tmatrix::RoundedCylinderGeometry<double> geom(half_height, radius);
                T_flat = tmatrix::compute_axisymmetric_tmatrix<double>(
                    geom, k, n_rel, lmax, Nint, use_ds, complex_plane, eps_z, conducting);
            }

            py::array_t<std::complex<double>> result({2, rmax, 2, rmax});
            auto buf = result.mutable_unchecked<4>();
            int idx = 0;
            for (int a1 = 0; a1 < 2; a1++)
                for (int r1 = 0; r1 < rmax; r1++)
                    for (int a2 = 0; a2 < 2; a2++)
                        for (int r2 = 0; r2 < rmax; r2++)
                            buf(a1, r1, a2, r2) = T_flat[idx++];

            return result;
        },
        py::arg("half_height"), py::arg("radius"),
        py::arg("k"), py::arg("n_rel"),
        py::arg("lmax"), py::arg("Nint") = 200,
        py::arg("use_ds") = true, py::arg("complex_plane") = false,
        py::arg("eps_z") = 0.95, py::arg("extended_precision") = false,
        py::arg("conducting") = false,
        R"pbdoc(
            Compute T-matrix for a rounded (oblate) cylinder using EBCM.

            Parameters
            ----------
            half_height : float
                Half-height of the cylinder.
            radius : float
                Radius including rounded part (must be > half_height).
            k : float
                Wavenumber in the medium.
            n_rel : complex
                Relative refractive index.
            lmax : int
                Maximum angular momentum order.
            Nint : int
                Number of quadrature points (default 200).
            use_ds : bool
                Use distributed sources (default True).
            complex_plane : bool
                Distribute sources in complex plane (default False).
            eps_z : float
                Source distribution extent parameter (default 0.95).
            extended_precision : bool
                Use __float128 precision internally (default False).
            conducting : bool
                If True, use PEC boundary conditions (default False).

            Returns
            -------
            T : ndarray, shape [2, rmax, 2, rmax]
        )pbdoc"
    );
}
