#include "interactions.hpp"
#include "vsh_translation.hpp"
#include "indices.hpp"
#include <math.h>
#include <Eigen/IterativeLinearSolvers>
#include <omp.h>

using std::complex;
using namespace std::complex_literals;

using Eigen::Vector3d;

ComplexVector matrix_vector_product(const Ref<const ComplexMatrix>& A, const Ref<const ComplexVector>& x) {
    int size = x.size();
    ComplexVector ret(size);
#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        ret(i) = A.row(i)*x;
    }

    return ret;
}

complex<double> dot_product(const Ref<const ComplexVector>& x, const Ref<const ComplexVector>& y) {
    int size = x.size();
    double real_part = 0;
    double imag_part = 0;
#pragma omp parallel for reduction(+:real_part,imag_part)
    for (int i = 0; i < size; i++) {
        real_part += x(i).real()*y(i).real() + x(i).imag()*y(i).imag();
        imag_part += x(i).real()*y(i).imag() - x(i).imag()*y(i).real();
    }

    return std::complex<double>(real_part, imag_part);
}

ComplexVector bicgstab(const Ref<const ComplexMatrix>& A, const Ref<const ComplexVector>& b,
        int maxiter, double tolerance) {

    int size = b.size();

    // step 1
    ComplexVector x_prev = b;

    double error = (A*x_prev - b).norm();
    if (error < tolerance)
        return x_prev;

    ComplexVector r_prev = b - A*x_prev;

    // step 2
    ComplexVector r_hat = r_prev;

    // step 3
    complex<double> rho_prev = 1;
    complex<double> alpha    = 1;
    complex<double> w_prev   = 1;

    // step 4
    ComplexVector v_prev = ComplexVector::Zero(size);
    ComplexVector p_prev = ComplexVector::Zero(size);

    // step 5
    int current_iteration = 1;

    while (true) {
        //complex<double> rho_i = r_hat.dot(r_prev);
        complex<double> rho_i = dot_product(r_hat, r_prev);
        complex<double> beta = (rho_i/rho_prev)*(alpha/w_prev);
        ComplexVector pi = r_prev + beta*(p_prev - w_prev*v_prev);
        ComplexVector vi = matrix_vector_product(A, pi);

        //alpha = rho_i/r_hat.dot(vi);
        alpha = rho_i/dot_product(r_hat, vi);
        ComplexVector h = x_prev + alpha*pi;

        ComplexVector s = r_prev - alpha*vi;
        ComplexVector t = matrix_vector_product(A, s);
        //complex<double> w_i = t.dot(s)/t.dot(t);
        complex<double> w_i = dot_product(t,s)/dot_product(t,t);
        ComplexVector xi = h + w_i*s;
        ComplexVector ri = s - w_i*t;

        error = ri.norm();
        if (error < tolerance || current_iteration > maxiter)
            return xi;

        x_prev   = xi;
        r_prev   = ri;
        rho_prev = rho_i;
        v_prev   = vi;
        p_prev   = pi;
        w_prev   = w_i;

        current_iteration += 1;
    }
}

ComplexVector solve_linear_system(const Ref<const ComplexMatrix>& agg_tmatrix,
        const Ref<const ComplexVector>& p_src, solver method) {

    ComplexMatrix interaction_matrix = agg_tmatrix;
    for (int i = 0; i < interaction_matrix.cols(); i++)
        interaction_matrix(i,i) += 1;
    
    switch (method) {
        case solver::bicgstab:
            //Eigen::BiCGSTAB<ComplexMatrix> solver;
            //solver.setTolerance(1e-5);
            //solver.compute(interaction_matrix);
            //return solver.solveWithGuess(p_src, p_src);
            return bicgstab(interaction_matrix, p_src);
    }
}

ComplexMatrix sphere_aggregate_tmatrix(const Ref<const position_t>& positions,
        const Ref<const ComplexMatrix>& mie, double k) {

    int Nparticles = positions.rows();
    int lmax = mie.cols()/2;
    int rmax = lmax_to_rmax(lmax);
    int size = 2*rmax*Nparticles;

    ComplexMatrix agg_tmatrix = ComplexMatrix::Zero(size, size);
    
    int N = Nparticles*(Nparticles-1)/2;
    Array ivals(N);
    Array jvals(N);
    int counter = 0;

    for (int i = 0; i < Nparticles; i++) {
        for (int j = i+1; j < Nparticles; j++) {
            ivals(counter) = i;
            jvals(counter) = j;
            counter += 1;
        }
    }

    if (Nparticles == 1)
        return agg_tmatrix;

    for (int n = 1; n < lmax+1; n++) {
        for (int m = -n; m < n+1; m++) {
            for (int v = 1; v < lmax+1; v++) {
                for (int u = -v; u < v+1; u++) {

                    auto fn = vsh_translation_lambda(m, n, u, v, vsh_mode::outgoing);

#pragma omp parallel for
                    for (int ij = 0; ij < N; ij++) {
                        int i = ivals(ij);
                        int j = jvals(ij);

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

    return agg_tmatrix;
}
