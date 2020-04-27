#include "decomposition.hpp"
#include "math.h"

using std::complex;
namespace py = pybind11;

complex<double> trapz(const Ref<const Array>& x, const Ref<const ComplexArray>& y) {
    int n = y.size();
    complex<double> sum = 0;

    for (int i=1; i < n; i++) {
        sum += (y(i) + y(i-1))*(x(i) - x(i-1))/2.0;
    }

    return sum;
}

complex<double> trapz_2d(const Ref<const Array>& x, const Ref<const Array>& y, const Ref<const ComplexMatrix>& f) {
    int n = y.size();
    ComplexArray z(n);

    for (int i=0; i < n; i++) {
        z(i) = trapz(x, f.col(i));
    }

    return trapz(y, z);
}

ComplexMatrix integrate_phase(const py::array_t<double> rhat, const Ref<const vec3>& origin, double k, int rmax,
        const Ref<const Array>& theta, const Ref<const Array>& phi, const py::array_t<std::complex<double>> p0) {

    auto rhat_p = rhat.unchecked<3>();
    auto p0_p = p0.unchecked<4>();

    const int Nx = theta.size();
    const int Ny = phi.size();

    ComplexMatrix p(2, rmax);

    ComplexMatrix exp_phase(Nx, Ny);
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double phase = -k*(rhat_p(0,i,j)*origin(0) + rhat_p(1,i,j)*origin(1) + rhat_p(2,i,j)*origin(2));
            exp_phase(i,j) = exp(complex<double>(0,phase));
        }
    }

    ComplexMatrix M(Nx, Ny);
    for (int i = 0; i < rmax; i++) {
        for (int a = 0; a < Nx; a++) {
            for (int b = 0; b < Ny; b++) {
                M(a,b) = p0_p(0,i,a,b)*exp_phase(a,b);
            }
        }
        p(0,i) = trapz_2d(theta, phi, M);

        for (int a = 0; a < Nx; a++) {
            for (int b = 0; b < Ny; b++) {
                M(a,b) = p0_p(1,i,a,b)*exp_phase(a,b);
            }
        }
        p(1,i) = trapz_2d(theta, phi, M);
    }

    return p;
}

using dtype = std::complex<double>;
using arr_in = pybind11::array_t<double>;
using arr_out = pybind11::array_t<dtype>;
arr_out grid_interpolate(const std::array<arr_in,2> grid, const arr_out data, const arr_in pts, dtype fill_value) {
    auto x_p = grid[0].unchecked<1>();
    auto y_p = grid[1].unchecked<1>();

    double x0_grid = x_p(0);
    double y0_grid = y_p(0);
    double dx_grid = x_p(1) - x_p(0);
    double dy_grid = y_p(1) - y_p(0);

    auto data_p = data.unchecked<3>();
    auto pts_p = pts.unchecked<2>();

    //auto pts_buf = pts.request();
    int Npts = pts_p.shape(0);
    int ndim_out = data.shape(2);
    auto result = arr_out({Npts, ndim_out});
    auto result_p = result.mutable_unchecked<2>();

    for (int i = 0; i < Npts; i++) {
        int x1 = std::floor((pts_p(i,0) - x0_grid) / dx_grid);
        int y1 = std::floor((pts_p(i,1) - y0_grid) / dy_grid);
        int x2 = x1 + 1;
        int y2 = y1 + 1;

        if (x1 < 0 || y1 < 0 || x2 > x_p.size() || y2 > y_p.size()) {
            for (int j = 0; j < ndim_out; j++) {
                result_p(i,j) = fill_value;
            }
        }
        else {
            double dx = pts_p(i,0) - x_p(x1);
            double dy = pts_p(i,1) - y_p(y1);

            for (int j = 0; j < ndim_out; j++) {
                result_p(i,j) = data_p(x1, y1, j)*(1 - dx/dx_grid)*(1 - dy/dy_grid)
                              + data_p(x2, y2, j)*(dx/dx_grid)*(dy/dy_grid)
                              + data_p(x1, y2, j)*(1 - dx/dx_grid)*(dy/dy_grid)
                              + data_p(x2, y1, j)*(dx/dx_grid)*(1 - dy/dy_grid);
            }
        }
    }

    return result;
}
