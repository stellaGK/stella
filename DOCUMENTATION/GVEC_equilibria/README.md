# stella–GVEC geometry tests

> **Location note.** This material was merged into the stella repository from the
> external GVEC working tree. Paths below are relative to the stella repo root.
>
> - `gvec_to_stella.py` is a **copy** kept here so the workflow is self-contained.
>   Its canonical home is the GVEC package (`gvec/python/gvec/scripts/gvec_to_stella.py`);
>   if you change one, change the other.
> - The GVEC/VMEC inputs and the comparison script live in this directory.

Two test cases for the new `geometry_option = 'gvec'` mode in stella.

---

## Test 1 — Circular tokamak (GVEC vs VMEC comparison)

The goal is to run both codes on an identical equilibrium and confirm that the
geometry arrays (bmag, grad\_alpha·grad\_alpha, etc.) agree to within ~1%.

### Files

| Path | Description |
|---|---|
| `circular_tok/gvec/parameter.ini` | GVEC input for circular tokamak (R0=3, a=1, q parabolic) |
| `circular_tok/vmec/input.tok_circ` | Matching VMEC input |

### Workflow

The automatic test `AUTOMATIC_TESTS/numerical_tests/test_2_geometry/test_2g_gvec_vs_vmec_geometry.py` runs stella with both geometries of this tokamak and of the W7-X equilibrium `wout_w7x_kjm.nc` (converted to GVEC with `pygvec convert-wout`), and compares the geometry arrays, the constants of the flux tube and a short linear ITG simulation. Its equilibrium files are made from the two inputs above by

```bash
pip install gvec==1.5.1 vmecpp==0.9.4 netCDF4
python AUTOMATIC_TESTS/numerical_tests/test_2_geometry/gvec_vs_vmec/make_gvec_vs_vmec_files.py <work_folder>
```

which solves the tokamak with VMEC++ and with GVEC, and writes the stella geometry files with `gvec_to_stella.py --a-minor <Aminor_p of the wout file>`, so that both geometries use the same reference length (GVEC's own `r_minor` is volume-based, VMEC's `Aminor_p` area-based). When stella is linked to the GVEC library, the test also reads the GVEC states directly (`gvec_state_file`; these runs are skipped otherwise), and it runs both paths with the magnetic field reversed (`gvec_vs_vmec/reverse_magnetic_field.py`), which must keep |B|, the metric, |φ|²(t), the growth rate and the heat flux, and reverse b·∇z, the drifts, the Clebsch factor and the real frequency, as for VMEC. The python interface of GVEC (`gvec_to_stella.py`) takes the toroidal flux from the parameter file (`phiedge`), while its Fortran interface (the direct path) takes the profiles stored in the state file, so a state must be used with the parameter file it was made with.

### Handedness

The logical coordinates of GVEC are right-handed, while the VMEC angles are left-handed (`signgs = -1`). The same equilibrium therefore has a toroidal flux and a rotational transform of opposite signs in the two codes: the VMEC input with `PHIEDGE > 0` and `iota > 0` corresponds to `PhiEdge < 0` and `iota < 0` in GVEC, as `pygvec convert-wout` gives and as `circular_tok/gvec/parameter.ini` sets. `gvec_to_stella.py` reverses the toroidal angle, which gives the VMEC angles for the cylindrical map of GVEC, and keeps the toroidal flux of GVEC as `psi_t`, as stella's VMEC geometry uses `psi_t = signgs * phi`. A frame of the wrong handedness gives the same `|B|`, metric and binormal drifts, but the opposite sign of the Clebsch factor, so the drive is reversed relative to the drifts.

### Expected result

With the resolutions of the two inputs, all geometry arrays agree to better than 0.5% at s = 0.25. Key scalar cross-checks at s = 0.25:

- `safety_factor_q` = 1.6 (q = 1.4 + 0.8 s)
- `shat` = 0.25 (= 2 s (dq/ds) / q)
- `L_reference` = 1.0 m (with `--a-minor 1.0`)
- `B_reference` = 1.0 T (= 2 |Phi_edge/(2π)| / a² = 2 × 0.5 / 1)

---

## Test 2 — Figure-8 (trefoil-knot axis stellarator)

This case uses GVEC's `which_hmap=10` (hmap\_knot) with a (2,3) torus-knot
axis — a figure-8-like topology that **VMEC cannot represent**.  The test
verifies that stella can read and use this exotic geometry without error.

### Files

| Path | Description |
|---|---|
| `figure8/gvec/parameter.ini` | GVEC input for trefoil-knot axis, nfp=3 |

### Workflow

```bash
# 1. Run GVEC
cd DOCUMENTATION/GVEC_equilibria/figure8/gvec
gvec parameter.ini
cd ../../../..

# 2. Post-process
python DOCUMENTATION/GVEC_equilibria/gvec_to_stella.py \
    DOCUMENTATION/GVEC_equilibria/figure8/gvec/GVEC_FIG8_State_00001.dat \
    --s 0.5 --nalpha 1 --nzgrid 192 --nfield_periods 3 \
    -o DOCUMENTATION/GVEC_equilibria/figure8/gvec_stella_fig8.nc

# 3. Self-consistency check only (no VMEC reference)
python DOCUMENTATION/GVEC_equilibria/compare_gvec_vmec.py \
    --gvec DOCUMENTATION/GVEC_equilibria/figure8/gvec_stella_fig8.nc
```

### Expected result

The self-consistency checks should all pass:
- `bmag > 0` everywhere
- `|∇ψt|² > 0` everywhere
- Cauchy-Schwarz inequality satisfied
- `q` finite and non-zero, `|shat|` reasonable

---

## Two ways to give stella a GVEC equilibrium

### 1. Direct (no pre-processing) — `gvec_state_file`

stella can read a GVEC state file itself, linking against the GVEC library, the
same way it reads a VMEC `wout`. No Python, no intermediate file, and the grid
can be changed freely without regenerating anything.

Build GVEC as a library with the gyrokinetic interface, then build stella
against it:

```bash
cmake <gvec-source> -DCOMPILE_GVEC_TO_GENE=ON -DCOMPILE_GVEC_AS_STATIC_LIB=ON
make -j
cmake <stella-source> -DSTELLA_ENABLE_GVEC=ON -DGVEC_ROOT=<gvec-build-dir>
make -j
```

or, with stella's Makefile build (the one used on Pitagora):

```bash
make USE_GVEC=on GVEC_ROOT=<gvec-build-dir> depend
make -j USE_GVEC=on GVEC_ROOT=<gvec-build-dir>
```

`GVEC_EXTRA_LIB` carries GVEC's BLAS/LAPACK dependency and defaults to
`-framework Accelerate` on macOS and `-llapack -lblas` elsewhere; override it in
a platform Makefile if your machine needs something else (e.g. MKL). Note that
`make clean` removes the CMake build directory as well, by design.

```fortran
&geometry_options
  geometry_option = 'gvec'
/
&geometry_gvec
  gvec_state_file     = 'MYCASE_State_final.dat'
  gvec_parameter_file = 'parameter_MYCASE_final.ini'
  torflux             = 0.5        ! required: which flux surface
/
```

`gvec_parameter_file` is needed for any equilibrium using a general (non-cylindrical)
frame mapping, because those h-maps read extra parameters of their own — the
axisNB frame reads `hmap_ncfile`, so that netCDF must also be reachable from the
run directory. Omitting this step is why GVEC's GENE interface is documented as
supporting only the cylindrical frame; with it, the general frames work.

`STELLA_ENABLE_GVEC` defaults to OFF, so stella still builds everywhere GVEC is
absent; selecting `gvec_state_file` in that case aborts with an explanatory
message rather than failing at link time.

### 2. Via `gvec_to_stella.py` — `gvec_filename`

The original route: pre-compute the geometry into a netCDF file. Still useful
when GVEC is not available on the machine running stella, or when an archivable
geometry file is wanted. Note that the file is tied to the stella grid — change
`nzed`, `ny`, `nfield_periods` or `z_grid_refinement_factor` and it must be
regenerated.

The two paths agree to better than 0.4% on the trefoil case (identical to
printed precision for `bmag`, `b.grad z` and `|grad x|^2`); the residual is the
difference between GVEC's SFL transform and pygvec's Newton solve for the PEST
angle.

## Stella input snippet

Once the geometry file is ready, add to your stella input file:

```fortran
&geometry_options
  geometry_option = 'gvec'
/
&geometry_gvec
  gvec_filename = 'path/to/gvec_stella_geometry.nc'
  z_grid_refinement_factor = 1      ! set to 4 if zed_equal_arc = .true.
  verbose = .true.
/
```
