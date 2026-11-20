# OpenMP Parallelisation Plan

## Background

stella uses hybrid MPI+OpenMP parallelisation (`MPI_THREAD_FUNNELED`).
MPI distributes work across two data layouts:

- **vmu_lo** — distributes `ivmu` (velocity-space index); used in `gyrokinetic_equation/`
- **kxkyz_lo** — distributes `ikxkyz` (k-space index); used in `dissipation/`, `fields/`

Each MPI rank owns a contiguous slice `llim_proc..ulim_proc` of its layout index.
OpenMP threads share that slice; each thread takes a subset of iterations.

### Conventions used throughout

```fortran
!$omp parallel default(none) &
!$omp firstprivate(layout_obj, scalar_grid_sizes) &
!$omp private(loop_indices, thread_local_work_arrays) &
!$omp shared(g, phi, pre_computed_read_only_arrays)
allocate (work_array(...))   ! per-thread allocation inside parallel region
!$omp do [collapse(N)]
do ivmu = ...
   ...
end do
!$omp end do
deallocate (work_array)
!$omp end parallel
```

Key rules:
- Always use `default(none)` — every variable must be explicit.
- Allocatable work arrays must be allocated **inside** the parallel region so each thread gets its own copy.
- `collapse(N)` requires perfectly nested loops (no statements between headers).
- Scalar assignments like `is = is_idx(...)` must move into the loop body when collapsing.

---

## Already done (earlier commits)

| File | What was parallelised |
|---|---|
| `calculations/calculations_add_explicit_terms.f90` | `add_explicit_term` and `add_explicit_term_ffs` |
| `gyrokinetic_equation/gk_mirror.f90` | `add_mirror_term`, `add_mirror_term_ffs`, FFS advance transform, `init_invert_mirror_operator` |
| `gyrokinetic_equation/gk_drive.f90` | FFS advance transform loop (`advance_wstar_explicit_ffs`) |
| `gyrokinetic_equation/gk_parallel_streaming.f90` | `get_dgdz`, `get_dgdz_centered` |
| `dissipation/collisions_fokkerplanck.f90` | `advance_collisions_fp_explicit`, both LU-solve loops in `advance_implicit_fp` |

---

## Tier A — High value, straightforward

### A1 · `gk_magnetic_drift.f90` — fluxtube gyro-average loops

**Runs:** every timestep  
**Subroutines:** `advance_wdrifty_explicit` (fluxtube branch), `advance_wdriftx_explicit` (fluxtube branch)  
**Pattern:** four serial `do ivmu` loops calling `gyro_average` / `gyro_average_j1`

```fortran
do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
    call gyro_average(dphidy, ivmu, g0k(:, :, :, :, ivmu))
end do
```

Each call writes only to `g0k(:,:,:,:,ivmu)` — a unique slice per thread.
`gyro_average` is stateless (no module-level writes).
The FFS paths in the same subroutines already have `!$omp`; these fluxtube paths are missing it.

**Fix:** `!$omp parallel do default(none) firstprivate(vmu_lo) private(ivmu) shared(dphidy, g0k)`

---

### A2 · `dissipation_hyper.f90` — apply hyper-dissipation

**Runs:** every timestep  
**Subroutine:** `advance_hyper_dissipation`

Two mutually exclusive branches (line ~232 and ~240):

- **`use_physical_ksqr` branch (line 232):** single `do ivmu` with an array expression — no private work arrays.
- **s-alpha branch (line 240):** `do ivmu / it / iz / iky` conditional divide — can collapse all four levels.

**Fix:** `!$omp parallel do` (outer loop only) with `shared(g, ...)`.

---

### A3 · `dissipation_hyper.f90` — `get_dgdvpa_fourth_order`

**Runs:** every timestep (inside `advance_hyper_vpa`)  
**Pattern:** `do ikxkyz` loop; `tmp(nvpa)` is a per-iteration scratch array.

```fortran
allocate (tmp(nvpa))
do ikxkyz = ...
    ...  ! fills tmp, then writes to gout(:, :, ikxkyz)
end do
deallocate (tmp)
```

`tmp` must become thread-private. Same pattern as `ghrs` in the FP implicit solver.

**Fix:** remove outer `allocate(tmp)`, wrap loop in `!$omp parallel`, allocate `tmp` inside, add `!$omp do`, deallocate before `!$omp end parallel`.

---

### A4 · `dissipation_hyper.f90` — `get_dgdz_fourth_order`

**Runs:** every timestep  
**Pattern:** `do ivmu` with inner `(iky, it)` loops; `gleft(2)` / `gright(2)` are fixed-size local scratch.

Identical structure to `get_dgdz` in `gk_parallel_streaming.f90` (already parallelised).

**Fix:** `!$omp parallel do default(none) collapse(2)` over `(ivmu, iky)` ... with `private(gleft, gright, ...)`.

---

## Tier B — Medium value, minor extra care needed

### B1 · `collisions_dougherty.f90` — `advance_collisions_dougherty_explicit`

**Runs:** every timestep  
**Pattern:** two `do ikxkyz` loops (one per `density_conservation` branch) calling `vpa_differential_operator` and `mu_differential_operator`; `mucoll(nmu)` is shared scratch.

Same structure as FP explicit loops already parallelised.  
`conserve_momentum` / `conserve_energy` are pure local-data routines — thread-safe.

**Fix:** make `mucoll` private (allocate inside parallel region); wrap each loop in `!$omp parallel ... !$omp do`.

---

### B2 · `gk_magnetic_drift.f90` — `init_wdrift` loops

**Runs:** once at initialisation  
**Subroutines:** `init_wdrift_without_neoclassical_terms` (line ~147), `init_wdrift_with_neoclassical_terms` (line ~256)

`wcvdrifty(nalpha, nzgrid)` and `wgbdrifty(nalpha, nzgrid)` are currently allocated once outside the loop and used as shared scratch. They must become per-thread private allocatables (same pattern as `energy` in `init_wstar`).

**Fix:** move `allocate(wcvdrifty/wgbdrifty/wcvdriftx/wgbdriftx)` inside `!$omp parallel`, add to `private()` clause.

---

### B3 · `collisions_dougherty.f90` — `init_vpadiff_matrix` / `init_mudiff_matrix`

**Runs:** once at initialisation  
**Pattern:** trivial `do ikxkyz` filling diagonal arrays `bb_vpa` / `bb_mu`. No private work arrays needed.

**Fix:** `!$omp parallel do default(none) firstprivate(...) private(ikxkyz, iky, ikx, iz, is) shared(bb_vpa/bb_mu, ...)`.

---

## Tier C — More complex, deferred

### C1 · `gk_parallel_nonlinearity.f90`

Three `do ivmu` loops (lines ~212, ~303, ~349) doing gyro-averages and Fourier transforms.
Thread-safe transforms are already ensured via `!$omp threadprivate` FFT buffers, but the loops
have several per-`ivmu` allocatables (`phi_gyro`, `dphidz`, `g0k`, `g0kxy`) that require
careful placement. Defer until Tier A/B are validated.

---

## Implementation order

1. [x] Tier A1 — `gk_magnetic_drift.f90` fluxtube gyro-average loops
2. [x] Tier A2 — `dissipation_hyper.f90` apply loops
3. [x] Tier A3 — `dissipation_hyper.f90` `get_dgdvpa_fourth_order`
4. [x] Tier A4 — `dissipation_hyper.f90` `get_dgdz_fourth_order`
5. [x] Tier B1 — `collisions_dougherty.f90` explicit loops
6. [x] Tier B2 — `gk_magnetic_drift.f90` init loops
7. [x] Tier B3 — `collisions_dougherty.f90` init matrix loops
8. [ ] Tier C1 — `gk_parallel_nonlinearity.f90` (deferred)
