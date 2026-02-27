# C++ EBCM T-Matrix Module — Notes for Future Work

## 1. Distributed Sources (DS) — Status: Disabled, Simple Fix Available

**Is this a regression?** Yes — the Fortran NFM-DS code uses DS by default for axisymmetric particles. The C++ code has DS disabled via `use_ds_internal = false` at `tmatrix_ebcm.cpp:248`.

**Root cause:** The DS coordinate vectors (`zRe`, `zIm`) are declared but never populated. `geom.distributed_sources()` is never called. When DS was tentatively enabled, `svwf_distributed()` accessed empty vectors → undefined behavior → garbage results.

**What's correct:** Everything *else* in the DS implementation matches the Fortran:
- Q-matrix dimension logic (NmaxL/NmaxC/Nrank) matches `Proces1.f90:298-305`
- Wavenumber arguments to `svwf_distributed()` match `mvnv_m()` in `MatrixQ.f90`
- SVWF computation matches `MN_DS` in `SVWF.f90:128-208`
- Mixed product assembly matches `matQ_m` for DS case

**Fix:** In `compute_axisymmetric_tmatrix()`, replace:
```cpp
bool use_ds_internal = false;
(void)complex_plane;
(void)eps_z;
std::vector<Real> zRe, zIm;
```
with:
```cpp
bool use_ds_internal = use_ds;
std::vector<Real> zRe, zIm;
if (use_ds_internal) {
    geom.distributed_sources(Nrank, complex_plane, eps_z, zRe, zIm);
}
```

**Why DS matters:** Localized sources converge slowly (or fail) for high-aspect-ratio particles (very elongated prolate or very flat oblate spheroids). DS distributes source points along the symmetry axis, improving conditioning. Without DS, extreme aspect ratios require higher Nint or may not converge at all.

**Validation after fix:** Compare C++ DS output against localized for moderate aspect ratios (2:1 spheroid should give same T-matrix either way). Then test extreme aspect ratios (10:1) where localized sources struggle.

---

## 2. Extended Precision (`__float128`) — Status: Works Internally, Inaccessible from High-Level API

**Infrastructure:** Complete and functional.
- `MIEPY_HAS_QUAD` properly detected (GCC on x86_64 Linux only)
- Eigen `NumTraits<__float128>` specialization complete
- All `tm_*` math functions have quadmath overloads
- All template instantiations present: `tmatrix_special.cpp`, `tmatrix_svwf.cpp`, `tmatrix_ebcm.cpp`
- pybind11 bindings route correctly via `#if MIEPY_HAS_QUAD`
- Downcast to double at output is explicit and correct

**What's broken:** All particle classes hardcode `extended_precision=False`:
- `spheroid.py:44`
- `cylinder.py` (both regular and rounded)
- No way for users to enable extended precision through the high-level API

**What's untested:**
- Zero tests for extended precision
- No test verifying quad gives better intermediates than double for challenging cases
- No test for platform fallback when `MIEPY_HAS_QUAD=0`
- Non-axisymmetric particles (ellipsoid, prisms) reference a `tmatrix_extended` Fortran binary that was never built — `extended_precision=True` would crash with `FileNotFoundError`

**Precision gain reality:** Output is always double (downcast), so final result precision is limited. The benefit is in intermediate calculations — fewer catastrophic cancellations for ill-conditioned problems. Difference between double and quad output is typically ~1e-16 (double rounding), but for challenging cases (high aspect ratio + large lmax), quad intermediates may prevent complete failure.

**To make accessible:**
- Add `extended_precision` parameter to particle constructors (or make it auto-detect based on condition number)
- Add test with a challenging case where double fails and quad succeeds
- Raise clear error for non-axisymmetric particles with `extended_precision=True`

---

## 3. Improvement Opportunities

### High Priority (high impact, low effort)

**A. Direct solve instead of explicit inverse** — `tmatrix_ebcm.cpp:218-220`
- Current: `Q31_inv = lu.inverse(); T_m = Q11 * Q31_inv;` (3x cubic ops)
- Better: `T_m = lu.solve(Q11);` (2x cubic ops, ~30% speedup, better numerics)
- 2-line change, no algorithmic deviation

**B. OpenMP parallelization of m-loop** — `tmatrix_ebcm.cpp:286`
- Each azimuthal mode m=0..lmax is independent (separate Q-matrix assembly + LU solve)
- Add `#pragma omp parallel for` with thread-local T arrays, merge at end
- Speedup proportional to cores (lmax+1 independent work items)
- Low effort, same algorithm

### Medium Priority (moderate impact or effort)

**C. SVWF allocation hoisting** — `tmatrix_svwf.cpp:20-21`
- `svwf_localized()` re-allocates `MV`/`NV` vectors (via `.assign()`) at every quadrature point
- For Nint=200, that's 800 heap allocations per Q-matrix per mode
- Fix: pre-allocate outside the quadrature loop in `assemble_Q_m()`, reuse storage
- ~5-15% speedup for the quadrature loop

**D. Error diagnostics for Q31 conditioning**
- No check if Q31 is singular or ill-conditioned → silent wrong results
- Add condition number estimate after LU; warn or fall back to SVD
- Prevents hard-to-debug failures for challenging parameter regimes

**E. Pybind11 code deduplication** — `tmatrix_pyb.cpp`
- 3 near-identical ~80-line blocks for spheroid/cylinder/rounded_cylinder
- Extract common reshape + dispatch logic into a template helper
- Code quality improvement, easier to add new geometries

### Low Priority (low impact or high effort)

**F. Quadrature loop OpenMP** — parallelizing the inner Nint loop
- Requires Eigen matrix reduction (no native `omp declare reduction` for Eigen)
- Need thread-local copies of A matrix + manual merge
- Medium effort for medium gain

**G. Adaptive quadrature** — auto-tune Nint for convergence
- Fortran has Mishchenko convergence criterion (in `MatrixQ.f90`)
- Would remove need for user to guess Nint
- Medium-high effort, deviates from current fixed-Nint interface

**H. Special function caching** — cache Bessel/Legendre evaluations
- Arguments vary continuously over the surface, so cache hit rate is low for localized sources
- Only valuable for DS mode where source coordinates are discrete
- Profile first before investing
