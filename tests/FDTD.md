# FDTD Validation Notes

## Running Tests

```bash
pixi run python tests/fdtd_validation.py --list                    # list test cases
pixi run python tests/fdtd_validation.py --run sphere_dielectric   # run one case
pixi run python tests/fdtd_validation.py --run all --resolution fast --plot
pixi run python tests/fdtd_validation.py --summary
```

## Resolution Presets

| Preset | Grid spacing | Time per test (3D) |
|--------|-------------|-------------------|
| `fast` | 12 nm | ~10-40 sec |
| `standard` | 8 nm | ~1-3 min |
| `high` | 5 nm | ~10-30 min |

## Current Results (fast, 12nm)

| Test | Scat L2 | Ext L2 |
|------|---------|--------|
| sphere_dielectric | 15.9% | 15.8% |
| spheroid_oblate | 7.1% | 7.2% |
| ellipsoid_aligned | 4.9% | 6.3% |
| cube_dielectric | 19.7% | 19.5% |
| hexprism_dielectric | 15.9% | 14.7% |

Errors are dominated by short-wavelength points (400-500nm) where sharp Mie resonances are under-resolved by FDTD. Mid-spectrum agreement is typically <5%.

## Followup: Higher Resolution

Standard and high resolution tests haven't been run for most cases yet. To get publishable-quality validation:

```bash
# Run all dielectric cases at standard (expect ~3-5% L2 for scat/ext)
pixi run python tests/fdtd_validation.py --run sphere_dielectric spheroid_oblate cylinder_aligned ellipsoid_aligned cube_dielectric hexprism_dielectric --resolution standard --plot

# Run gold cases at standard
pixi run python tests/fdtd_validation.py --run sphere_gold spheroid_gold --resolution standard --plot

# High resolution for sphere baseline (expect <2% L2)
pixi run python tests/fdtd_validation.py --run sphere_dielectric --resolution high --plot
```

High resolution is expensive: the 5nm grid spacing means 8x more voxels than fast (12nm), so expect ~8-30x longer runtimes. Each 3D test at high resolution takes 10-30 minutes.

Cached results persist in `tests/fdtd_cache/` across runs. Use `--force` to re-run.

## Followup: High Aspect Ratio Tests

The `spheroid_prolate` (5:1) and `spheroid_extreme` (10:1) tests stress both FDTD and T-matrix:

**FDTD side**: The thin axis (80nm for 5:1, 60nm for 10:1) needs enough grid points to resolve. At 12nm fast resolution, the 60nm axis gets only 5 pixels — not enough. These cases need at minimum standard (8nm, ~8 pixels) and ideally high (5nm, ~12 pixels) resolution.

```bash
# Prolate 5:1 — should work at standard resolution
pixi run python tests/fdtd_validation.py --run spheroid_prolate --resolution standard --plot

# Extreme 10:1 — needs high resolution for the 60nm thin axis
pixi run python tests/fdtd_validation.py --run spheroid_extreme --resolution high --plot
```

**T-matrix side**: The extreme spheroid uses `extended_precision=True` in miepy and lmax=5. The EBCM solver may still warn about ill-conditioned Q matrices at high m-values. If results look wrong, try increasing lmax in `_miepy_spheroid_extreme()`.

**Expected errors**: Prolate 5:1 should be <10% at standard resolution. The 10:1 extreme case may show 10-20% errors even at high resolution — both methods struggle with extreme aspect ratios, so larger discrepancies are acceptable.

## Followup: Rotated Particle Tests

`cylinder_tilted` and `ellipsoid_rotated` haven't been run yet. These test the orientation/rotation handling in both Meep (axis/e1/e2/e3 parameters) and miepy (quaternion-based T-matrix rotation):

```bash
pixi run python tests/fdtd_validation.py --run cylinder_tilted ellipsoid_rotated --resolution fast --plot
```

If these show larger errors than their aligned counterparts, the rotation mapping between miepy quaternions and Meep's e1/e2/e3 vectors may need debugging. The mapping is: `R = quaternion.as_rotation_matrix(q)`, then `e1=R[:,0]`, `e2=R[:,1]`, `e3=R[:,2]`.
