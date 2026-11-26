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

```bash
# 1. Run GVEC
cd DOCUMENTATION/GVEC_equilibria/circular_tok/gvec
gvec parameter.ini
cd ../../../..

# 2. Post-process GVEC → stella geometry file
python DOCUMENTATION/GVEC_equilibria/gvec_to_stella.py \
    DOCUMENTATION/GVEC_equilibria/circular_tok/gvec/GVEC_TOK_CIRC_State_00001.dat \
    --s 0.5 --nalpha 1 --nzgrid 128 \
    -o DOCUMENTATION/GVEC_equilibria/circular_tok/gvec_stella_tok_circ.nc

# 3. Run VMEC (using your local VMEC installation)
cd DOCUMENTATION/GVEC_equilibria/circular_tok/vmec
xvmec input.tok_circ
cd ../../../..

# 4. Post-process VMEC wout → same stella geometry format
python DOCUMENTATION/GVEC_equilibria/gvec_to_stella.py \
    DOCUMENTATION/GVEC_equilibria/circular_tok/vmec/wout_tok_circ.nc \
    --s 0.5 --nalpha 1 --nzgrid 128 \
    -o DOCUMENTATION/GVEC_equilibria/circular_tok/vmec_stella_tok_circ.nc
# (this requires a separate vmec_to_stella.py script or the VMEC branch of stella)

# Alternatively: let stella itself read both and compare its diagnostic output.

# 5. Compare
python DOCUMENTATION/GVEC_equilibria/compare_gvec_vmec.py \
    --gvec DOCUMENTATION/GVEC_equilibria/circular_tok/gvec_stella_tok_circ.nc \
    --vmec DOCUMENTATION/GVEC_equilibria/circular_tok/vmec_stella_tok_circ.nc \
    --tol 1e-2
```

### Expected result

All geometry quantities should agree within ~1% (tolerance is set by the
difference in resolution and the minor differences between GVEC's B-spline
basis and VMEC's Fourier series).  Key scalar cross-checks:

- `safety_factor_q` ≈ 1.63 at s=0.5 (q0=1.4 + 0.8*0.5)
- `shat` ≈ 0.49 (= -rho * diota/drho / iota at rho=sqrt(0.5))
- `L_reference` ≈ 1.0 m
- `B_reference` ≈ 2.0 T (= 2|Phi_edge/(2π)|/a² = 2*0.5/1 = 1.0 T)

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
- `q > 0`, `|shat|` reasonable

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
