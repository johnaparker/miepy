#include "decomposition.hpp"

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
