#include "flux.hpp"
#include "indices.hpp"
#include <cmath>

std::array<Array,3> particle_cross_sections(const Ref<const ComplexVector>& p_scat, const Ref<const ComplexVector>& p_inc,
        const Ref<const ComplexVector>& p_src, double k) {

    int rmax = p_scat.size()/2;
    int lmax = rmax_to_lmax(rmax);

    Array Cext = Array::Zero(2*lmax);
    Array Cabs = Array::Zero(2*lmax);

    double factor = 4*PI/pow(k,2);

    for (int n = 1; n < lmax + 1; n++) {
        for (int m = -n; m < n+1; m++) {
            for (int a = 0; a < 2; a++) {
                int r = a*rmax + pow(n,2) + n + m - 1;
                int idx = a*lmax + n - 1;
                Cabs[idx] += factor*(real(conj(p_inc[r])*p_scat[r]) - norm(p_scat[r]));
                Cext[idx] += factor*real(conj(p_src[r])*p_scat[r]);
            }
        }
    }

    Array Cscat = Cext - Cabs;

    return std::array<Array,3>{Cscat, Cabs, Cext};
}

std::array<Array,3> cluster_cross_sections(const Ref<const ComplexVector>& p_cluster,
        const Ref<const ComplexVector>& p_src, double k) {

    int rmax = p_cluster.size()/2;
    int lmax = rmax_to_lmax(rmax);

    Array Cscat = Array::Zero(2*lmax);
    Array Cext  = Array::Zero(2*lmax);

    double factor = 4*PI/pow(k,2);

    for (int n = 1; n < lmax + 1; n++) {
        for (int m = -n; m < n+1; m++) {
            for (int a = 0; a < 2; a++) {
                int r = a*rmax + pow(n,2) + n + m - 1;
                int idx = a*lmax + n - 1;
                Cscat[idx] += factor*norm(p_cluster[r]);
                Cext[idx]  += factor*real(conj(p_src[r])*p_cluster[r]);
            }
        }
    }

    Array Cabs = Cext - Cscat;

    return std::array<Array,3>{Cscat, Cabs, Cext};
}
