#include "interactions.hpp"
#include <math.h>

using std::complex;
using namespace std::complex_literals;

using Eigen::Vector3d;

int lmax_to_rmax(int lmax) {
    return lmax*(lmax + 2);
}

ComplexMatrix sphere_aggregate_tmatrix(const Ref<const position_t>& positions,
        const Ref<const ComplexMatrix>& mie, double k) {

    int Nparticles = positions.rows();
    int lmax = mie.cols()/2;
    int rmax = lmax_to_rmax(lmax);
    int size = 2*rmax*Nparticles;

    ComplexMatrix agg_tmatrix = ComplexMatrix::Zero(size, size);
    
    if (Nparticles == 1)
        return agg_tmatrix;

    for (int n = 1; n < lmax+1; n++) {
        for (int m = -n; m < n+1; m++) {
            for (int v = 1; v < lmax+1; v++) {
                for (int u = -v; u < v+1; u++) {

                    auto fn = vsh_translation_lambda(m, n, u, v, vsh_mode::outgoing);

                    for (int i = 0; i < Nparticles; i++) {
                        for (int j = i+1; j < Nparticles; j++) {
                            Vector3d dji = positions.row(i) - positions.row(j);

                            double rad = dji.norm();
                            double theta = acos(dji(2)/rad);
                            double phi = atan2(dji(1), dji(0));
                            auto transfer = fn(rad, theta, phi, k);

                            for (int a = 0; a < 2; a++) {
                                for (int b = 0; b < 2; b++) {
                                    complex<double> val = transfer[(a+b)%2];
                                    int idx = i*(2*rmax) + a*(rmax) + n*(n+2) - n + m - 1;
                                    int idy = j*(2*rmax) + b*(rmax) + v*(v+2) - v + u - 1;

                                    agg_tmatrix(idx, idy) = val*mie(j, b*lmax + v-1);
                                    agg_tmatrix(idy, idx) = pow(-1, n+v+a+b)*val*mie(i, a*lmax + n-1);
                                }
                            }
                        }
                    }
                }
            }
        }
    } 

    return agg_tmatrix;
}
