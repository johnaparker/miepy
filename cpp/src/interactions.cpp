#define NOMINMAX
#include "interactions.hpp"
#include "indices.hpp"
#include "vsh_translation.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <omp.h>

using std::complex;
using namespace std::complex_literals;

using Eigen::Vector3d;

ComplexVector matrix_vector_product(const Ref<const ComplexMatrix> &A,
                                    const Ref<const ComplexVector> &x) {
  int size = x.size();
  ComplexVector ret(size);
#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    ret(i) = (A.row(i) * x).value();
  }

  return ret;
}

complex<double> dot_product(const Ref<const ComplexVector> &x,
                            const Ref<const ComplexVector> &y) {
  int size = x.size();
  double real_part = 0;
  double imag_part = 0;
#pragma omp parallel for reduction(+ : real_part, imag_part)
  for (int i = 0; i < size; i++) {
    real_part += x(i).real() * y(i).real() + x(i).imag() * y(i).imag();
    imag_part += x(i).real() * y(i).imag() - x(i).imag() * y(i).real();
  }

  return std::complex<double>(real_part, imag_part);
}

ComplexVector bicgstab(const Ref<const ComplexMatrix> &A,
                       const Ref<const ComplexVector> &b, int maxiter,
                       double tolerance) {

  int size = b.size();

  // step 1
  ComplexVector x_prev = b;
  ComplexVector r_prev = b - matrix_vector_product(A, x_prev);

  double error = (r_prev).norm();
  if (error < tolerance)
    return x_prev;

  // step 2
  ComplexVector r_hat = r_prev;

  // step 3
  complex<double> rho_prev = 1;
  complex<double> alpha = 1;
  complex<double> w_prev = 1;

  // step 4
  ComplexVector v_prev = ComplexVector::Zero(size);
  ComplexVector p_prev = ComplexVector::Zero(size);

  // step 5
  int current_iteration = 1;

  while (true) {
    // complex<double> rho_i = r_hat.dot(r_prev);
    complex<double> rho_i = dot_product(r_hat, r_prev);
    complex<double> beta = (rho_i / rho_prev) * (alpha / w_prev);
    ComplexVector pi = r_prev + beta * (p_prev - w_prev * v_prev);
    ComplexVector vi = matrix_vector_product(A, pi);

    // alpha = rho_i/r_hat.dot(vi);
    alpha = rho_i / dot_product(r_hat, vi);
    ComplexVector h = x_prev + alpha * pi;

    ComplexVector s = r_prev - alpha * vi;
    ComplexVector t = matrix_vector_product(A, s);
    // complex<double> w_i = t.dot(s)/t.dot(t);
    complex<double> w_i = dot_product(t, s) / dot_product(t, t);
    ComplexVector xi = h + w_i * s;
    ComplexVector ri = s - w_i * t;

    error = ri.norm();
    if (error < tolerance || current_iteration > maxiter)
      return xi;

    x_prev = xi;
    r_prev = ri;
    rho_prev = rho_i;
    v_prev = vi;
    p_prev = pi;
    w_prev = w_i;

    current_iteration += 1;
  }
}

bicgstab_result bicgstab_profiled(const Ref<const ComplexMatrix> &A,
                       const Ref<const ComplexVector> &b, int maxiter,
                       double tolerance) {

  int size = b.size();

  // step 1
  ComplexVector x_prev = b;
  ComplexVector r_prev = b - matrix_vector_product(A, x_prev);

  double error = (r_prev).norm();
  if (error < tolerance)
    return {x_prev, 0, error};

  // step 2
  ComplexVector r_hat = r_prev;

  // step 3
  complex<double> rho_prev = 1;
  complex<double> alpha = 1;
  complex<double> w_prev = 1;

  // step 4
  ComplexVector v_prev = ComplexVector::Zero(size);
  ComplexVector p_prev = ComplexVector::Zero(size);

  // step 5
  int current_iteration = 1;

  while (true) {
    complex<double> rho_i = dot_product(r_hat, r_prev);
    complex<double> beta = (rho_i / rho_prev) * (alpha / w_prev);
    ComplexVector pi = r_prev + beta * (p_prev - w_prev * v_prev);
    ComplexVector vi = matrix_vector_product(A, pi);

    alpha = rho_i / dot_product(r_hat, vi);
    ComplexVector h = x_prev + alpha * pi;

    ComplexVector s = r_prev - alpha * vi;
    ComplexVector t = matrix_vector_product(A, s);
    complex<double> w_i = dot_product(t, s) / dot_product(t, t);
    ComplexVector xi = h + w_i * s;
    ComplexVector ri = s - w_i * t;

    error = ri.norm();
    if (error < tolerance || current_iteration > maxiter)
      return {xi, current_iteration, error};

    x_prev = xi;
    r_prev = ri;
    rho_prev = rho_i;
    v_prev = vi;
    p_prev = pi;
    w_prev = w_i;

    current_iteration += 1;
  }
}

ComplexVector apply_block_preconditioner(
    const Ref<const ComplexMatrix> &M_inv_blocks,
    const Ref<const ComplexVector> &x, int block_size) {

  int size = x.size();
  int N = size / block_size;
  ComplexVector result(size);

#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    int offset = i * block_size;
    result.segment(offset, block_size) =
        M_inv_blocks.block(offset, 0, block_size, block_size) *
        x.segment(offset, block_size);
  }

  return result;
}

ComplexVector bicgstab_preconditioned(const Ref<const ComplexMatrix> &A,
                                      const Ref<const ComplexVector> &b,
                                      const Ref<const ComplexMatrix> &M_inv_blocks,
                                      int block_size, int maxiter,
                                      double tolerance) {

  int size = b.size();

  ComplexVector x_prev = b;
  ComplexVector r_prev = b - matrix_vector_product(A, x_prev);

  double error = r_prev.norm();
  if (error < tolerance)
    return x_prev;

  ComplexVector r_hat = r_prev;

  complex<double> rho_prev = 1;
  complex<double> alpha = 1;
  complex<double> w_prev = 1;

  ComplexVector v_prev = ComplexVector::Zero(size);
  ComplexVector p_prev = ComplexVector::Zero(size);

  int current_iteration = 1;

  while (true) {
    complex<double> rho_i = dot_product(r_hat, r_prev);
    complex<double> beta = (rho_i / rho_prev) * (alpha / w_prev);
    ComplexVector pi = r_prev + beta * (p_prev - w_prev * v_prev);

    // Right preconditioning: p_hat = M^{-1} * p, then v = A * p_hat
    ComplexVector p_hat = apply_block_preconditioner(M_inv_blocks, pi, block_size);
    ComplexVector vi = matrix_vector_product(A, p_hat);

    alpha = rho_i / dot_product(r_hat, vi);
    ComplexVector h = x_prev + alpha * p_hat;

    ComplexVector s = r_prev - alpha * vi;

    // Right preconditioning: s_hat = M^{-1} * s, then t = A * s_hat
    ComplexVector s_hat = apply_block_preconditioner(M_inv_blocks, s, block_size);
    ComplexVector t = matrix_vector_product(A, s_hat);

    complex<double> w_i = dot_product(t, s) / dot_product(t, t);
    ComplexVector xi = h + w_i * s_hat;
    ComplexVector ri = s - w_i * t;

    error = ri.norm();
    if (error < tolerance || current_iteration > maxiter)
      return xi;

    x_prev = xi;
    r_prev = ri;
    rho_prev = rho_i;
    v_prev = vi;
    p_prev = pi;
    w_prev = w_i;

    current_iteration += 1;
  }
}

bicgstab_result bicgstab_preconditioned_profiled(
    const Ref<const ComplexMatrix> &A, const Ref<const ComplexVector> &b,
    const Ref<const ComplexMatrix> &M_inv_blocks, int block_size, int maxiter,
    double tolerance) {

  int size = b.size();

  ComplexVector x_prev = b;
  ComplexVector r_prev = b - matrix_vector_product(A, x_prev);

  double error = r_prev.norm();
  if (error < tolerance)
    return {x_prev, 0, error};

  ComplexVector r_hat = r_prev;

  complex<double> rho_prev = 1;
  complex<double> alpha = 1;
  complex<double> w_prev = 1;

  ComplexVector v_prev = ComplexVector::Zero(size);
  ComplexVector p_prev = ComplexVector::Zero(size);

  int current_iteration = 1;

  while (true) {
    complex<double> rho_i = dot_product(r_hat, r_prev);
    complex<double> beta = (rho_i / rho_prev) * (alpha / w_prev);
    ComplexVector pi = r_prev + beta * (p_prev - w_prev * v_prev);

    ComplexVector p_hat = apply_block_preconditioner(M_inv_blocks, pi, block_size);
    ComplexVector vi = matrix_vector_product(A, p_hat);

    alpha = rho_i / dot_product(r_hat, vi);
    ComplexVector h = x_prev + alpha * p_hat;

    ComplexVector s = r_prev - alpha * vi;

    ComplexVector s_hat = apply_block_preconditioner(M_inv_blocks, s, block_size);
    ComplexVector t = matrix_vector_product(A, s_hat);

    complex<double> w_i = dot_product(t, s) / dot_product(t, t);
    ComplexVector xi = h + w_i * s_hat;
    ComplexVector ri = s - w_i * t;

    error = ri.norm();
    if (error < tolerance || current_iteration > maxiter)
      return {xi, current_iteration, error};

    x_prev = xi;
    r_prev = ri;
    rho_prev = rho_i;
    v_prev = vi;
    p_prev = pi;
    w_prev = w_i;

    current_iteration += 1;
  }
}

ComplexVector solve_linear_system_preconditioned(
    const Ref<const ComplexMatrix> &agg_tmatrix,
    const Ref<const ComplexVector> &p_src,
    const Ref<const ComplexMatrix> &M_inv_blocks, int block_size,
    solver method) {

  ComplexMatrix interaction_matrix = agg_tmatrix;
  for (int i = 0; i < interaction_matrix.cols(); i++)
    interaction_matrix(i, i) += 1;

  switch (method) {
  case solver::exact:
  default:
  case solver::bicgstab:
    return bicgstab_preconditioned(interaction_matrix, p_src, M_inv_blocks,
                                   block_size);
  }
}

ComplexVector solve_linear_system(const Ref<const ComplexMatrix> &agg_tmatrix,
                                  const Ref<const ComplexVector> &p_src,
                                  solver method) {

  ComplexMatrix interaction_matrix = agg_tmatrix;
  for (int i = 0; i < interaction_matrix.cols(); i++)
    interaction_matrix(i, i) += 1;

  switch (method) {
  case solver::exact: // TODO: implment
  default:
  case solver::bicgstab:
    // Eigen::BiCGSTAB<ComplexMatrix> solver;
    // solver.setTolerance(1e-5);
    // solver.compute(interaction_matrix);
    // return solver.solveWithGuess(p_src, p_src);
    return bicgstab(interaction_matrix, p_src);
  }
}

ComplexMatrix sphere_aggregate_tmatrix(const Ref<const position_t> &positions,
                                       const Ref<const ComplexMatrix> &mie,
                                       double k) {

  int lmax = mie.cols() / 2;
  int rmax = lmax_to_rmax(lmax);
  int Nparticles = positions.rows();
  int size = 2 * rmax * Nparticles;

  ComplexMatrix agg_tmatrix = ComplexMatrix::Zero(size, size);

  if (Nparticles == 1)
    return agg_tmatrix;

  int N = Nparticles * (Nparticles - 1) / 2;
  Array ivals(N);
  Array jvals(N);
  int counter = 0;

  for (int i = 0; i < Nparticles; i++) {
    for (int j = i + 1; j < Nparticles; j++) {
      ivals(counter) = i;
      jvals(counter) = j;
      counter += 1;
    }
  }

  auto vsh_precompute = create_vsh_cache_map(lmax);

#pragma omp parallel for
  for (int ij = 0; ij < N; ij++) {
    int i = ivals(ij);
    int j = jvals(ij);

    Vector3d dji = positions.row(i) - positions.row(j);

    double rad = dji.norm();
    double theta = acos(dji(2) / rad);
    double phi = atan2(dji(1), dji(0));

    vsh_translation_insert_pair(agg_tmatrix, mie, i, j, rad, theta, phi, k,
                                vsh_precompute);
  }

  return agg_tmatrix;
}

ComplexMatrix particle_aggregate_tmatrix(const Ref<const position_t> &positions,
                                         const tmatrix_t &tmatrix, double k) {

  int rmax = tmatrix.dimensions()[1] / 2;
  int lmax = rmax_to_lmax(rmax);

  int Nparticles = positions.rows();
  int size = 2 * rmax * Nparticles;

  ComplexMatrix agg_tmatrix = ComplexMatrix::Zero(size, size);

  if (Nparticles == 1)
    return agg_tmatrix;

  int N = Nparticles * (Nparticles - 1) / 2;
  Array ivals(N);
  Array jvals(N);
  int counter = 0;

  for (int i = 0; i < Nparticles; i++) {
    for (int j = i + 1; j < Nparticles; j++) {
      ivals(counter) = i;
      jvals(counter) = j;
      counter += 1;
    }
  }

  auto vsh_precompute = create_vsh_cache_map(lmax);

#pragma omp parallel for
  for (int ij = 0; ij < N; ij++) {
    int i = ivals(ij);
    int j = jvals(ij);

    Vector3d dji = positions.row(i) - positions.row(j);

    double rad = dji.norm();
    double theta = acos(dji(2) / rad);
    double phi = atan2(dji(1), dji(0));

    vsh_translation_insert_pair(agg_tmatrix, tmatrix, i, j, rad, theta, phi, k,
                                vsh_precompute);
  }

  return agg_tmatrix;
}

ComplexMatrix reflection_matrix_nia(const Ref<const position_t> &positions,
                                    const Ref<const ComplexMatrix> &mie,
                                    double k, complex<double> reflection,
                                    double z) {

  int lmax = mie.cols() / 2;
  int rmax = lmax_to_rmax(lmax);
  int Nparticles = positions.rows();
  int size = 2 * rmax * Nparticles;

  ComplexMatrix R_matrix = ComplexMatrix::Zero(size, size);

  for (int i = 0; i < Nparticles; i++) {
    Vector3d pi = positions.row(i);
    pi(2) -= 2 * (pi(2) - z);
    for (int j = 0; j < Nparticles; j++) {
      Vector3d pj = positions.row(j);
      Vector3d dji = pi - pj;

      double rad = dji.norm();
      double theta = acos(dji(2) / rad);
      double phi = atan2(dji(1), dji(0));

      for (int n = 1; n < lmax + 1; n++) {
        for (int m = -n; m < n + 1; m++) {
          for (int v = 1; v < lmax + 1; v++) {
            for (int u = -v; u < v + 1; u++) {
              auto transfer = vsh_translation(m, n, u, v, rad, theta, phi, k,
                                              vsh_mode::outgoing);
              for (int a = 0; a < 2; a++) {
                for (int b = 0; b < 2; b++) {
                  complex<double> factor = -reflection * pow(-1, m + n + a + 1);
                  complex<double> val = transfer[(a + b) % 2];
                  int idx =
                      i * (2 * rmax) + a * (rmax) + n * (n + 2) - n + m - 1;
                  int idy =
                      j * (2 * rmax) + b * (rmax) + v * (v + 2) - v + u - 1;
                  R_matrix(idx, idy) = factor * val * mie(j, b * lmax + v - 1);
                }
              }
            }
          }
        }
      }
    }
  }

  return R_matrix;
}
