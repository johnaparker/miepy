#ifndef GUARD_vsh_translation_h
#define GUARD_vsh_translation_h

#include <complex>
#include <functional>
#include <tuple>

enum class vsh_mode {
    outgoing,
    ingoing,
    incident,
    interior
};

std::function<std::complex<double>(int, double, bool)> get_zn(vsh_mode mode);

std::tuple<std::complex<double>, std::complex<double>> vsh_translation(
        int m, int n, int u, int v, double rad, double theta,
        double phi, double k, vsh_mode mode);

double test3();

#endif
