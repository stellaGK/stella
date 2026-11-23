#!/usr/bin/env python3
"""compare_stella_geometry.py
=============================
Validate stella's GVEC geometry branch.

Two modes:

``--check FILE``
    Run physical-invariant checks on a single ``<run>.geometry`` file:
    |B| > 0, a positive-definite perpendicular metric, and the Cauchy-Schwarz
    inequality. These hold for any valid equilibrium, so this works for the
    figure-8 and trefoil cases where there is no VMEC reference.

``--compare A B``
    Compare two ``<run>.geometry`` files element by element. The intended use
    is a cross-code benchmark: run the SAME equilibrium through stella twice,
    once with ``geometry_option = 'vmec'`` (reading the wout directly) and once
    with ``geometry_option = 'gvec'``, then compare.

Why compare stella's own output rather than the input files?
------------------------------------------------------------
Everything downstream of the reader - the (alpha, psi_t) -> (x, y) transform,
the boundary-condition factors, the z-grid - is shared between the VMEC and
GVEC branches. Comparing ``.geometry`` files therefore isolates the difference
that actually matters: how each code represents the same equilibrium. It also
needs no extra converter, because stella reads a VMEC wout natively.

Usage
-----
::

    python compare_stella_geometry.py --check  trefoil.geometry
    python compare_stella_geometry.py --compare w7x_vmec.geometry w7x_gvec.geometry --tol 1e-2
"""

import argparse
import sys

import numpy as np

# Column layout of stella's <run>.geometry file, in order. See
# write_geometric_coefficients() in STELLA_CODE/geometry/geometry.f90.
COLUMNS = [
    "alpha", "zed", "zeta", "bmag", "b_dot_gradz",
    "grady_dot_grady", "gradx_dot_grady", "gradx_dot_gradx",
    "BxgradB_dot_grady", "Bxkappa_dot_grady", "BxgradB_dot_gradx",
    "bmag_psi0",
]

# Quantities compared in --compare mode. 'zed'/'zeta'/'alpha' are grid
# coordinates and are checked separately for alignment.
COMPARED = [
    "bmag", "b_dot_gradz", "grady_dot_grady", "gradx_dot_grady",
    "gradx_dot_gradx", "BxgradB_dot_grady", "Bxkappa_dot_grady",
    "BxgradB_dot_gradx",
]


def read_geometry(path):
    """Parse a stella <run>.geometry file into (scalars, arrays)."""
    header, rows = [], []
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                header.append(stripped)
                continue
            rows.append([float(v) for v in stripped.split()])

    if not rows:
        sys.exit(f"ERROR: no data rows found in {path}")

    data = np.array(rows)
    if data.shape[1] != len(COLUMNS):
        sys.exit(
            f"ERROR: {path} has {data.shape[1]} columns, expected {len(COLUMNS)}.\n"
            "       The .geometry format may have changed; update COLUMNS."
        )
    arrays = {name: data[:, i] for i, name in enumerate(COLUMNS)}

    # The second header line holds the scalar flux-surface quantities.
    scalars = {}
    scalar_names = ["rhoc", "qinp", "shat", "rhotor", "aref", "bref",
                    "dxdpsi", "dydalpha", "exb_nonlin", "flux_fac", "one_over_Grho"]
    for line in header:
        parts = line.lstrip("#").split()
        if len(parts) == len(scalar_names):
            try:
                values = [float(p) for p in parts]
            except ValueError:
                continue
            scalars = dict(zip(scalar_names, values))
            break
    return scalars, arrays


def check_invariants(path, rtol):
    """Physical invariants that must hold for any valid equilibrium."""
    scalars, g = read_geometry(path)
    failures = []

    print(f"\nChecking {path}")
    if scalars:
        print(f"  q = {scalars['qinp']:.5f}   shat = {scalars['shat']:.5f}   "
              f"a_ref = {scalars['aref']:.5f} m   B_ref = {scalars['bref']:.5f} T")
    print(f"  {len(g['bmag'])} grid points, "
          f"zeta in [{g['zeta'].min():.4f}, {g['zeta'].max():.4f}]")

    def require(condition, message):
        status = "ok  " if condition else "FAIL"
        print(f"    [{status}] {message}")
        if not condition:
            failures.append(message)

    require(np.all(g["bmag"] > 0), "|B| > 0 everywhere")
    require(np.all(np.isfinite(np.column_stack([g[c] for c in COLUMNS]))),
            "all geometry entries are finite")
    require(np.all(g["gradx_dot_gradx"] > 0), "|grad x|^2 > 0 everywhere")
    require(np.all(g["grady_dot_grady"] > 0), "|grad y|^2 > 0 everywhere")

    lhs = g["gradx_dot_grady"] ** 2
    rhs = g["gradx_dot_gradx"] * g["grady_dot_grady"]
    require(np.all(lhs <= rhs * (1.0 + rtol)),
            "Cauchy-Schwarz: (grad x . grad y)^2 <= |grad x|^2 |grad y|^2")

    # b . grad z must not change sign: the field line is followed monotonically.
    require(np.all(g["b_dot_gradz"] > 0) or np.all(g["b_dot_gradz"] < 0),
            "b . grad z does not change sign")

    if failures:
        print(f"\n  {len(failures)} invariant(s) FAILED")
    else:
        print("\n  all invariants satisfied")
    return not failures


def compare(path_a, path_b, tol):
    """Element-by-element comparison of two .geometry files."""
    scalars_a, a = read_geometry(path_a)
    scalars_b, b = read_geometry(path_b)

    print(f"\nComparing\n  A = {path_a}\n  B = {path_b}")

    if len(a["bmag"]) != len(b["bmag"]):
        sys.exit(
            f"ERROR: different number of grid points "
            f"({len(a['bmag'])} vs {len(b['bmag'])}).\n"
            "       Use the same nzed / nfield_periods for both runs."
        )

    if scalars_a and scalars_b:
        print("\n  Flux-surface scalars:")
        print(f"    {'quantity':<14}{'A':>14}{'B':>14}{'rel diff':>12}")
        for name in ["qinp", "shat", "aref", "bref"]:
            va, vb = scalars_a[name], scalars_b[name]
            rel = abs(va - vb) / max(abs(va), abs(vb), 1e-30)
            flag = "" if rel <= tol else "   <-- exceeds tol"
            print(f"    {name:<14}{va:>14.6e}{vb:>14.6e}{rel:>12.2e}{flag}")

    # The two runs must sit on the same z grid for a pointwise comparison.
    dz = np.max(np.abs(a["zed"] - b["zed"]))
    print(f"\n  max |zed_A - zed_B| = {dz:.3e}")
    if dz > 1e-8:
        print("    WARNING: the z grids differ; the pointwise comparison below "
              "is not meaningful.")

    print("\n  Geometry arrays (relative L2 and max-norm differences):")
    print(f"    {'quantity':<22}{'rel L2':>12}{'rel max':>12}{'':>6}")
    worst = 0.0
    failures = []
    for name in COMPARED:
        va, vb = a[name], b[name]
        scale = max(np.linalg.norm(va), np.linalg.norm(vb), 1e-30)
        rel_l2 = np.linalg.norm(va - vb) / scale
        denom = np.maximum(np.abs(va), np.abs(vb))
        # Compare against the array scale, not pointwise, so near-zero
        # crossings do not dominate the metric.
        rel_max = np.max(np.abs(va - vb)) / max(np.max(denom), 1e-30)
        worst = max(worst, rel_l2)
        flag = "ok" if rel_l2 <= tol else "FAIL"
        if rel_l2 > tol:
            failures.append(name)
        print(f"    {name:<22}{rel_l2:>12.3e}{rel_max:>12.3e}{flag:>6}")

    print(f"\n  worst relative L2 difference = {worst:.3e}   (tolerance {tol:.1e})")
    if failures:
        print(f"  {len(failures)} quantity(ies) exceed the tolerance: {', '.join(failures)}")
    else:
        print("  all compared quantities agree within tolerance")
    return not failures


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", metavar="FILE",
                        help="run invariant checks on a single .geometry file")
    parser.add_argument("--compare", nargs=2, metavar=("A", "B"),
                        help="compare two .geometry files")
    parser.add_argument("--tol", type=float, default=1e-2,
                        help="relative tolerance (default: 1e-2)")
    args = parser.parse_args()

    if not args.check and not args.compare:
        parser.error("give --check FILE or --compare A B")

    ok = True
    if args.check:
        ok &= check_invariants(args.check, args.tol)
    if args.compare:
        ok &= compare(args.compare[0], args.compare[1], args.tol)

    print()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
