# OpenMP Parallelisation Plan — Field Equations Solver

## Background

The `STELLA_CODE/field_equations/` directory contains **no OpenMP directives** at all.
All loops over MPI-local index ranges are currently serial.

The field equations are solved using two MPI data layouts:

- **vmu_lo** layout — `g(naky, nakx, -nzgrid:nzgrid, ntubes, ivmu)` with `ivmu` MPI-distributed.
  Used by `*_vmulo` subroutine variants.
- **kxkyz_lo** layout — `g(nvpa, nmu, ikxkyz)` with `ikxkyz = (iky, ikx, iz, it, is)` MPI-distributed.
  Used by `*_kxkyzlo` subroutine variants.

In the `kxkyz_lo` layout a single index `ikxkyz` encodes the full tuple
`(iky, ikx, iz, it, is)`. This has a critical consequence for thread-safety
(see Category 2 below).

---

## Category 1 — Race-free, straightforward to implement

These loops write to **unique array elements per iteration** — no atomics needed.

---

### C1-A · `field_equations_electromagnetic.f90` — vmulo mu/vpa multiply loops

**Runs:** every timestep (electromagnetic simulations)
**Subroutine:** `advance_fields_electromagnetic_vmulo`

```fortran
! Line 184 — multiply phi_gyro by mu factor
do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
   imu = imu_idx(vmu_lo, ivmu)
   phi_gyro(:, :, :, :, ivmu) = phi_gyro(:, :, :, :, ivmu) * mu(imu)
end do

! Line 219 — multiply phi_gyro by vpa factor
do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
   iv = iv_idx(vmu_lo, ivmu)
   phi_gyro(:, :, :, :, ivmu) = phi_gyro(:, :, :, :, ivmu) * vpa(iv)
end do
```

Each iteration writes to a unique `ivmu` slice of `phi_gyro`. No private work arrays needed.

**Fix:** `!$omp parallel do default(none) firstprivate(vmu_lo) private(ivmu, imu/iv) shared(phi_gyro, mu/vpa)`

---

### C1-B · `field_equations_electromagnetic.f90` — `calculate_phi_and_bpar`

**Runs:** every timestep (electromagnetic + bpar simulations)
**Subroutine:** `calculate_phi_and_bpar` (lines 441–452)

```fortran
do it = 1, ntubes
   do iz = -nzgrid, nzgrid
      do ikx = 1, nakx
         do iky = 1, naky
            antot1 = phi(iky,ikx,iz,it)
            antot3 = bpar(iky,ikx,iz,it)
            phi(iky,ikx,iz,it)  = inv11(iky,ikx,iz)*antot1 + inv13(iky,ikx,iz)*antot3
            bpar(iky,ikx,iz,it) = inv31(iky,ikx,iz)*antot1 + inv33(iky,ikx,iz)*antot3
         end do
      end do
   end do
end do
```

Each `(it, iz, ikx, iky)` is read and written exactly once. Four perfectly nested loops,
no branching between headers.

**Fix:** `!$omp parallel do default(none) collapse(4) private(it, iz, ikx, iky, antot1, antot3) shared(phi, bpar, inv11, ...)`

---

### C1-C · `field_equations_electromagnetic.f90` — init matrix inverse loops

**Runs:** once at initialisation
**Subroutine:** `init_field_equations_electromagnetic` (lines 758–781)

Triple nested `do iz / do ikx / do iky` computing
`denominator_fields_inv{11,13,31,33}` element-wise from `denominator_fields`.
Each `(iz, ikx, iky)` is independent.

**Fix:** `!$omp parallel do default(none) collapse(3) private(iz, ikx, iky, denom_tmp) shared(denominator_fields*, ...)`

---

### C1-D · `field_equations_electromagnetic.f90` — init kxkyz accumulation loops (×3)

**Runs:** once at initialisation
**Subroutine:** `init_field_equations_electromagnetic` (lines 688–702, 713–728, 734–749)

Three `do ikxkyz` loops filling `apar_denom`, `denominator_fields33`, and
`denominator_fields13`. Each contains `if (it /= 1) cycle`, so only the
`it == 1` slice of `ikxkyz` is active. After the cycle, each remaining
iteration writes to a unique `(iky, ikx, iz)` location.

`g0(nvpa, nmu)` is a local temporary — needs to be private (allocate inside `!$omp parallel`).

**Fix:** wrap each loop in `!$omp parallel ... !$omp do` with `g0` private.

---

### C1-E · `field_equations_fluxtube.f90` — init denominator loop

**Runs:** once at initialisation
**Subroutine:** `init_field_equations_fluxtube` (lines 689–709)

Same pattern as C1-D: `do ikxkyz` with `if (it /= 1) cycle`, writing to
`denominator_fields(iky, ikx, iz)`. After the cycle, unique grid point per iteration.

`g0(nvpa, nmu)` local temporary — make private.

**Fix:** same pattern as C1-D.

---

### C1-F · `field_equations_collisions.f90` — `get_fields_by_spec`

**Runs:** every timestep (Fokker-Planck collision operator)
**Subroutine:** `get_fields_by_spec` (lines 84–94)

```fortran
do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
   ...
   call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
   g0 = g0 * wgt
   call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
end do
```

The output array `fld` is 5D including the species index `is`. Each `ikxkyz`
encodes a unique `(iky, ikx, iz, it, is)` tuple — no two iterations write the
same element. `g0(nvpa, nmu)` must be private.

**Fix:** `!$omp parallel ... !$omp do` with `g0` private (allocate per thread).

---

### C1-G · `field_equations_collisions.f90` — `get_fields_by_spec_idx`

**Runs:** every timestep (Fokker-Planck collision operator)
**Subroutine:** `get_fields_by_spec_idx` (lines 177–191)

Same output pattern as C1-F (`fld(..., is)` is unique per `ikxkyz`), but with
an inline Bessel function computation into `g0(:, imu)` inside a `do imu` loop.
`g0(nvpa, nmu)` must be private.

**Fix:** same as C1-F.

---

## Category 2 — Species-reduction race (requires atomics or restructuring)

These are the primary timestep field-solve kernels and the hottest loops in the
field solver. They accumulate species contributions into `phi`/`apar`/`bpar`:

```fortran
do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
   iky = iky_idx(...); ikx = ikx_idx(...); iz = iz_idx(...); it = it_idx(...)
   is  = is_idx(...)
   ...
   phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp  ! ← RACE
end do
```

Two `ikxkyz` values that differ only in species `is` map to the **same**
`(iky, ikx, iz, it)` grid point. Both accumulate into `phi`. Parallelizing the
`ikxkyz` loop without synchronisation produces wrong results.

Affected subroutines:
- `advance_fields_fluxtube_using_field_equations_kxkyzlo` — `phi +=` (line 342)
- `advance_fields_electromagnetic_kxkyzlo` — `phi +=`, `bpar +=` (lines 313, 319), `apar +=` (line 363)
- `advance_apar` — `apar +=` (line 581)

### Options

**Option A — `!$omp atomic` on the update line**

```fortran
!$omp atomic update
phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
```

Requires OpenMP 5.1 for complex-valued atomics. Not all Fortran compilers
support this yet (Intel ifx ≥ 2023.x, GFortran ≥ 13 with `-fopenmp`).

**Option B — Parallelize the spatial `(iky, ikx)` loops instead**

After the `kxkyz_lo` loop (which stays serial), the division step
`phi = phi / denominator_fields` is a full-grid loop that can be parallelized
with `collapse(4)` over `(it, iz, ikx, iky)`. This captures some of the
available parallelism without race conditions.

**Option C — Restructure to explicit species loop (larger refactor)**

Split the single `do ikxkyz` into an outer `do is` loop and an inner spatial
loop. The inner loop is then race-free. This is the cleanest solution but
requires significant restructuring of the kxkyz_lo layout logic.

**Recommendation:** implement Option B now (easy win on the division step), defer
Options A/C until compiler support is confirmed or the refactor is budgeted.

---

## Implementation order

| # | Item | File | Lines | Runs | Effort |
|---|---|---|---|---|---|
| 1 | C1-A vmulo mu/vpa multiply | `field_equations_electromagnetic.f90` | 184, 219 | every step | trivial |
| 2 | C1-B `calculate_phi_and_bpar` | `field_equations_electromagnetic.f90` | 441–452 | every step | trivial |
| 3 | C1-F `get_fields_by_spec` | `field_equations_collisions.f90` | 84–94 | every step | easy |
| 4 | C1-G `get_fields_by_spec_idx` | `field_equations_collisions.f90` | 177–191 | every step | easy |
| 5 | C1-C init matrix inverse | `field_equations_electromagnetic.f90` | 758–781 | once | easy |
| 6 | C1-D init kxkyz loops (×3) | `field_equations_electromagnetic.f90` | 688–749 | once | medium |
| 7 | C1-E init denominator | `field_equations_fluxtube.f90` | 689–709 | once | medium |
| 8 | C2 Option B division step | `field_equations_fluxtube.f90`, `field_equations_electromagnetic.f90` | various | every step | easy |
| 9 | C2 Options A or C (species sum) | multiple | various | every step | complex |

Items 1–4 are the highest priority: called every timestep, race-free, and require
minimal changes. Items 5–7 are straightforward init-time wins. Items 8–9 address
the core field-solve bottleneck but require more care.
