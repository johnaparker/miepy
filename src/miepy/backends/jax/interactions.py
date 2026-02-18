"""JAX implementation of aggregate T-matrix assembly and linear system solver.

The T-matrix assembly uses a two-phase design:
  Phase A: Precompute all mode tuples, Gaunt coefficients, and scatter indices (NumPy, cached per lmax)
  Phase B: JIT-compiled JAX core that vectorizes over particle pairs and modes (runs on GPU)
"""

import functools
from math import factorial

import numpy as np

import miepy


def _precompute_mode_data(lmax):
    """Precompute all static index arrays and Gaunt coefficients for given lmax.

    Returns a dict of NumPy arrays ready for JAX conversion.
    Called once per lmax value and cached by the caller.
    """
    from miepy.backends.jax.special import gaunt_batch

    rmax = lmax * (lmax + 2)

    # Step 1: Enumerate valid mode tuples (n, m, v, u) with v <= n,
    # skipping the duplicate case n==v and u < -m
    n_list, m_list, v_list, u_list = [], [], [], []

    for n in range(1, lmax + 1):
        for m in range(-n, n + 1):
            for v in range(1, n + 1):
                for u in range(-v, v + 1):
                    if n == v and u < -m:
                        continue
                    n_list.append(n)
                    m_list.append(m)
                    v_list.append(v)
                    u_list.append(u)

    M = len(n_list)
    n_arr = np.array(n_list, dtype=np.int32)
    m_arr = np.array(m_list, dtype=np.int32)
    v_arr = np.array(v_list, dtype=np.int32)
    u_arr = np.array(u_list, dtype=np.int32)
    neg_m_arr = -m_arr

    # Step 2: Prefactors and Gaunt coefficient data
    factor = np.zeros(M, dtype=np.float64)
    mu_sum = (u_arr + neg_m_arr).astype(np.float64)

    gaunt_data = []
    for t in range(M):
        n, v = int(n_arr[t]), int(v_arr[t])
        u, neg_m = int(u_arr[t]), int(neg_m_arr[t])

        factor[t] = 0.5 * (-1)**neg_m * np.sqrt(
            (2*v + 1) * (2*n + 1) * factorial(v - u) * factorial(n - neg_m)
            / (v*(v+1) * n*(n+1) * factorial(v + u) * factorial(n + neg_m))
        )

        gaunt_a, gaunt_b = gaunt_batch(neg_m, n, u, v)
        qmax_A = min(n, v, (n + v - abs(neg_m + u)) // 2)
        qmax_B = min(n, v, (n + v + 1 - abs(neg_m + u)) // 2)
        gaunt_data.append((gaunt_a, gaunt_b, qmax_A, qmax_B))

    # Step 3: Padded Gaunt coefficient arrays with masks
    Q = max(max(gd[2] + 1 for gd in gaunt_data),
            max(gd[3] + 1 for gd in gaunt_data))

    A_coeffs = np.zeros((M, Q), dtype=np.complex128)
    B_coeffs = np.zeros((M, Q), dtype=np.complex128)
    A_mask = np.zeros((M, Q), dtype=np.float64)
    B_mask = np.zeros((M, Q), dtype=np.float64)
    A_pnm_idx = np.zeros((M, Q), dtype=np.int32)
    A_zn_idx = np.zeros((M, Q), dtype=np.int32)
    B_pnm_idx = np.zeros((M, Q), dtype=np.int32)
    B_zn_idx = np.zeros((M, Q), dtype=np.int32)

    for t in range(M):
        n, v = int(n_arr[t]), int(v_arr[t])
        mu = int(mu_sum[t])
        gaunt_a, gaunt_b, qmax_A, qmax_B = gaunt_data[t]

        for q in range(qmax_A + 1):
            p = n + v - 2*q
            A_coeffs[t, q] = gaunt_a[q] * (1j)**p * (n*(n+1) + v*(v+1) - p*(p+1))
            A_mask[t, q] = 1.0
            A_pnm_idx[t, q] = p*(p + 2) - p + mu
            A_zn_idx[t, q] = p

        for q in range(1, qmax_B + 1):
            p = n + v - 2*q
            B_coeffs[t, q] = gaunt_b[q] * (1j)**(p+1) * np.sqrt(
                ((p+1)**2 - (n-v)**2) * ((n+v+1)**2 - (p+1)**2)
            )
            B_mask[t, q] = 1.0
            B_pnm_idx[t, q] = (p+1)*(p+3) - (p+1) + mu
            B_zn_idx[t, q] = p + 1

    # Step 4: Scatter index arrays for all 4 entry types x 4 (a,b) combos
    scatter_row_local = []
    scatter_col_local = []
    scatter_row_is_i = []
    scatter_col_is_i = []
    scatter_sign = []
    scatter_mode_idx = []
    scatter_ab_parity = []
    scatter_mie_is_i = []
    scatter_mie_pol = []
    scatter_mie_order = []

    def _add(row_l, col_l, ri, ci, sgn, tidx, abp, mi, mp, mo):
        scatter_row_local.append(row_l)
        scatter_col_local.append(col_l)
        scatter_row_is_i.append(ri)
        scatter_col_is_i.append(ci)
        scatter_sign.append(sgn)
        scatter_mode_idx.append(tidx)
        scatter_ab_parity.append(abp)
        scatter_mie_is_i.append(mi)
        scatter_mie_pol.append(mp)
        scatter_mie_order.append(mo)

    for t in range(M):
        n, m, v, u = int(n_arr[t]), int(m_arr[t]), int(v_arr[t]), int(u_arr[t])
        has_34 = (n != v) or (u != -m)

        for a in range(2):
            for b in range(2):
                abp = (a + b) % 2

                # Entry 1: i <- j
                _add(a*rmax + n*(n+2) - n + m - 1,
                     b*rmax + v*(v+2) - v + u - 1,
                     True, False, 1.0, t, abp, False, b, v - 1)

                # Entry 2: j <- i (sign flip)
                _add(a*rmax + n*(n+2) - n + m - 1,
                     b*rmax + v*(v+2) - v + u - 1,
                     False, True, float((-1)**(n+v+a+b)), t, abp, True, b, v - 1)

                if has_34:
                    # Entry 3: negated m,u, i <- j
                    _add(b*rmax + v*(v+2) - v - u - 1,
                         a*rmax + n*(n+2) - n - m - 1,
                         True, False, float((-1)**(m+u-a-b)), t, abp, False, a, n - 1)

                    # Entry 4: negated m,u, j <- i
                    _add(b*rmax + v*(v+2) - v - u - 1,
                         a*rmax + n*(n+2) - n - m - 1,
                         False, True, float((-1)**(m+u+n+v)), t, abp, True, a, n - 1)

    return {
        'rmax': rmax,
        'factor': factor,
        'mu_sum': mu_sum,
        'A_coeffs': A_coeffs,
        'B_coeffs': B_coeffs,
        'A_mask': A_mask,
        'B_mask': B_mask,
        'A_pnm_idx': A_pnm_idx,
        'A_zn_idx': A_zn_idx,
        'B_pnm_idx': B_pnm_idx,
        'B_zn_idx': B_zn_idx,
        'scatter_row_local': np.array(scatter_row_local, dtype=np.int32),
        'scatter_col_local': np.array(scatter_col_local, dtype=np.int32),
        'scatter_row_is_i': np.array(scatter_row_is_i),
        'scatter_col_is_i': np.array(scatter_col_is_i),
        'scatter_sign': np.array(scatter_sign, dtype=np.float64),
        'scatter_mode_idx': np.array(scatter_mode_idx, dtype=np.int32),
        'scatter_ab_parity': np.array(scatter_ab_parity, dtype=np.int32),
        'scatter_mie_is_i': np.array(scatter_mie_is_i),
        'scatter_mie_pol': np.array(scatter_mie_pol, dtype=np.int32),
        'scatter_mie_order': np.array(scatter_mie_order, dtype=np.int32),
    }


@functools.lru_cache(maxsize=16)
def _get_sphere_agg_tmatrix_jit(lmax):
    """Get a JIT-compiled sphere_aggregate_tmatrix function for given lmax.

    Precomputed arrays become XLA constants baked into the compiled program.
    One JIT compilation per (lmax, N) pair (N from input shapes).
    """
    import jax
    import jax.numpy as jnp

    from miepy.backends.jax.special import (
        associated_legendre_recursion,
        spherical_hn_recursion,
    )

    precomp = _precompute_mode_data(lmax)
    rmax = precomp['rmax']
    p_max = 2 * lmax + 1
    block_size = 2 * rmax

    # Convert precomputed arrays to JAX (become XLA constants in JIT closure)
    factor = jnp.asarray(precomp['factor'])
    mu_sum = jnp.asarray(precomp['mu_sum'])
    A_coeffs = jnp.asarray(precomp['A_coeffs'])
    B_coeffs = jnp.asarray(precomp['B_coeffs'])
    A_mask = jnp.asarray(precomp['A_mask'])
    B_mask = jnp.asarray(precomp['B_mask'])
    A_pnm_idx = jnp.asarray(precomp['A_pnm_idx'])
    A_zn_idx = jnp.asarray(precomp['A_zn_idx'])
    B_pnm_idx = jnp.asarray(precomp['B_pnm_idx'])
    B_zn_idx = jnp.asarray(precomp['B_zn_idx'])
    s_row_local = jnp.asarray(precomp['scatter_row_local'])
    s_col_local = jnp.asarray(precomp['scatter_col_local'])
    s_row_is_i = jnp.asarray(precomp['scatter_row_is_i'])
    s_col_is_i = jnp.asarray(precomp['scatter_col_is_i'])
    s_sign = jnp.asarray(precomp['scatter_sign'])
    s_mode_idx = jnp.asarray(precomp['scatter_mode_idx'])
    s_ab_parity = jnp.asarray(precomp['scatter_ab_parity'])
    s_mie_is_i = jnp.asarray(precomp['scatter_mie_is_i'])
    s_mie_pol = jnp.asarray(precomp['scatter_mie_pol'])
    s_mie_order = jnp.asarray(precomp['scatter_mie_order'])

    @jax.jit
    def core(positions, mie, k, i_vals, j_vals):
        N = positions.shape[0]
        size = N * block_size

        # Step 1: Pair geometry [P = N*(N-1)/2 pairs]
        dji = positions[i_vals] - positions[j_vals]         # [P, 3]
        rad = jnp.linalg.norm(dji, axis=1)                  # [P]
        cos_theta = dji[:, 2] / rad                          # [P]
        phi = jnp.arctan2(dji[:, 1], dji[:, 0])             # [P]

        # Step 2: Special functions via vmap (Python loops unrolled at trace time)
        zn_all = jax.vmap(
            lambda r: spherical_hn_recursion(p_max, k * r)
        )(rad)                                                # [P, p_max+1]
        Pnm_all = jax.vmap(
            lambda ct: associated_legendre_recursion(p_max, ct)
        )(cos_theta)                                          # [P, pnm_size]

        # Step 3: Vectorized translation computation
        exp_phi = jnp.exp(1j * mu_sum * phi[:, None])         # [P, M]

        # A translations: sum_q A_coeff[q] * Pnm[idx_A[q]] * zn[p_A[q]]
        Pnm_A = Pnm_all[:, A_pnm_idx]                        # [P, M, Q]
        zn_A = zn_all[:, A_zn_idx]                             # [P, M, Q]
        sum_A = jnp.sum(A_coeffs * Pnm_A * zn_A * A_mask, axis=-1)  # [P, M]
        A_trans = factor * exp_phi * sum_A                      # [P, M]

        # B translations: -factor * exp_phi * sum_q B_coeff[q] * Pnm[idx_B[q]] * zn[p_B[q]]
        Pnm_B = Pnm_all[:, B_pnm_idx]                        # [P, M, Q]
        zn_B = zn_all[:, B_zn_idx]                             # [P, M, Q]
        sum_B = jnp.sum(B_coeffs * Pnm_B * zn_B * B_mask, axis=-1)  # [P, M]
        B_trans = -factor * exp_phi * sum_B                     # [P, M]

        transfer = jnp.stack([A_trans, B_trans], axis=-1)       # [P, M, 2]

        # Step 4: Scatter into flat output matrix
        # Determine particle index for row/col/mie of each (pair, entry)
        row_particle = jnp.where(s_row_is_i, i_vals[:, None], j_vals[:, None])  # [P, E]
        col_particle = jnp.where(s_col_is_i, i_vals[:, None], j_vals[:, None])  # [P, E]
        mie_particle = jnp.where(s_mie_is_i, i_vals[:, None], j_vals[:, None])  # [P, E]

        rows = row_particle * block_size + s_row_local         # [P, E]
        cols = col_particle * block_size + s_col_local         # [P, E]

        transfer_vals = transfer[:, s_mode_idx, s_ab_parity]   # [P, E]
        mie_vals = mie[mie_particle, s_mie_pol, s_mie_order]  # [P, E]
        vals = s_sign * transfer_vals * mie_vals                # [P, E]

        # All (row, col) indices are globally unique -- .set() is safe.
        # Scatter real and imaginary parts separately as float64 to avoid
        # XLA's slow complex128 scatter on GPU (no native complex atomics).
        flat_idx = rows.ravel() * size + cols.ravel()
        flat_vals = vals.ravel()
        agg_real = jnp.zeros(size * size, dtype=jnp.float64).at[flat_idx].set(flat_vals.real)
        agg_imag = jnp.zeros(size * size, dtype=jnp.float64).at[flat_idx].set(flat_vals.imag)
        agg_flat = agg_real + 1j * agg_imag

        return agg_flat.reshape(N, 2, rmax, N, 2, rmax)

    return core


def sphere_aggregate_tmatrix_jax(positions, mie, k):
    """Compute the sphere aggregate T-matrix using fully vectorized JAX.

    Arguments:
        positions[N,3]   particle positions
        mie[N,2,lmax]    Mie scattering coefficients
        k                medium wavenumber

    Returns:
        T[N,2,rmax,N,2,rmax] aggregate T-matrix (JAX array)
    """
    import jax.numpy as jnp

    positions = jnp.asarray(positions, dtype=jnp.float64)
    mie = jnp.asarray(mie, dtype=jnp.complex128)
    N = positions.shape[0]
    lmax = mie.shape[-1]
    rmax = lmax * (lmax + 2)

    if N == 1:
        return jnp.zeros((1, 2, rmax, 1, 2, rmax), dtype=jnp.complex128)

    jit_fn = _get_sphere_agg_tmatrix_jit(lmax)
    i_vals, j_vals = np.triu_indices(N, k=1)
    i_vals = jnp.asarray(i_vals)
    j_vals = jnp.asarray(j_vals)
    return jit_fn(positions, mie, k, i_vals, j_vals)


def _make_jit_solvers():
    """Create JIT-compiled solver functions (traced once per matrix shape)."""
    import jax
    import jax.numpy as jnp
    import jax.scipy.sparse.linalg as jsla

    @jax.jit
    def _solve_exact(A, b):
        return jnp.linalg.solve(A, b)

    @jax.jit
    def _solve_bicgstab(A, b):
        matvec = lambda v: A @ v
        x, info = jsla.bicgstab(matvec, b, tol=1e-5, maxiter=1000)
        return x

    return _solve_exact, _solve_bicgstab

_jit_solvers = None

def _get_jit_solvers():
    global _jit_solvers
    if _jit_solvers is None:
        _jit_solvers = _make_jit_solvers()
    return _jit_solvers


def solve_linear_system_jax(agg_tmatrix, p_src, method):
    """Solve the linear system p_inc = (I + T_agg)^{-1} * p_src using JAX.

    Arguments:
        agg_tmatrix[N,2,rmax,N,2,rmax]   aggregate T-matrix
        p_src[N,2,rmax]                   source coefficients
        method                            solver method (miepy.solver.exact or .bicgstab)

    Returns:
        p_inc[N,2,rmax]
    """
    import jax.numpy as jnp

    Nparticles = agg_tmatrix.shape[0]
    rmax = p_src.shape[-1]
    size = Nparticles * 2 * rmax

    # Convert to JAX arrays and flatten
    T_2d = jnp.asarray(agg_tmatrix.reshape(size, size))
    b = jnp.asarray(p_src.reshape(-1))

    # Build (I + T) matrix
    A = jnp.eye(size, dtype=jnp.complex128) + T_2d

    solve_exact, solve_bicgstab = _get_jit_solvers()
    if method == miepy.solver.exact:
        x = solve_exact(A, b)
    else:
        x = solve_bicgstab(A, b)

    return np.asarray(x).reshape(Nparticles, 2, rmax)
