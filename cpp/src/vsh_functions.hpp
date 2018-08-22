#ifndef GUARD_vsh_functions_h
#define GUARD_vsh_functions_h

#include <functional>
#include <complex>
#include "vec.hpp"

using E_type = Eigen::Matrix<std::complex<double>, 3, Eigen::Dynamic, Eigen::RowMajor>;

enum class vsh_mode {
    outgoing,
    ingoing,
    incident,
    interior
};

std::function<std::complex<double>(int, double, bool)> get_zn(vsh_mode mode);

std::complex<double> Emn(int m, int n);

cvec3 vsh_electric(int n, int m, vsh_mode mode, double rad,
        double theta, double phi, double k);

cvec3 vsh_magnetic(int n, int m, vsh_mode mode, double rad,
        double theta, double phi, double k);

E_type expand_E_cluster(const Ref<const position_t>& pos, const Ref<const ComplexMatrix>& p,
        vsh_mode mode, const Ref<const Array>& x, const Ref<const Array>& y,
        const Ref<const Array>& z, double k);

#endif
