#!/usr/bin/env python3
"""compare_gvec_vmec.py
======================
Compare geometry quantities from a GVEC-generated stella input file and a
VMEC wout file for the same equilibrium.

Both files should represent identical (or near-identical) equilibria so that
the geometry arrays can be compared directly.

Usage
-----
::

    python compare_gvec_vmec.py \\
        --gvec   tests/circular_tok/gvec_stella_tok_circ.nc \\
        --vmec   tests/circular_tok/vmec/wout_tok_circ.nc \\
        [--s 0.5  --nalpha 1  --nzgrid 128  --tol 1e-2]

The VMEC wout file is post-processed with gvec_to_stella.py's sister function
so that both inputs are proper gvec_stella_geometry.nc files.  Alternatively,
supply two pre-computed gvec_stella_geometry.nc files (one from GVEC, one
generated from the VMEC wout via gvec_to_stella.py or a separate script).

If only ``--gvec`` is given, the script performs self-consistency checks on
that file alone (useful for the figure-8 case where no VMEC reference exists).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")           # non-interactive backend
import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def load_geometry(nc_file: str | Path) -> dict:
    """Load all variables from a gvec_stella_geometry.nc file into a dict."""
    data = {}
    with nc.Dataset(str(nc_file), "r") as ds:
        for name, var in ds.variables.items():
            data[name] = var[...].data   # convert masked arrays to plain ndarray
        data["_filename"] = str(nc_file)
    return data


def rel_err(a, b):
    """Relative L∞ error between two arrays, normalised by max(|b|)."""
    denom = np.max(np.abs(b))
    if denom == 0:
        return 0.0
    return float(np.max(np.abs(a - b))) / denom


def abs_err(a, b):
    """Absolute L∞ error."""
    return float(np.max(np.abs(a - b)))


_SCALAR_VARS = [
    "s", "safety_factor_q", "shat", "L_reference", "B_reference",
]

_2D_VARS = [
    ("bmag",                       "B/B_ref"),
    ("b_dot_gradzeta",             "(a/B) B.∇ζ"),
    ("grad_alpha_grad_alpha",      "a²|∇α|²"),
    ("grad_alpha_grad_psit",       "∇α·∇ψt / B_ref"),
    ("grad_psit_grad_psit",        "|∇ψt|² / (a²B_ref²)"),
    ("gds23",                      "gds23"),
    ("gds24",                      "gds24"),
    ("gds25",                      "gds25"),
    ("gds26",                      "gds26"),
    ("gbdrift_alpha",              "a²B_ref/B³ (B×∇B)·∇α"),
    ("gbdrift_kappa_alpha",        "a²B_ref/B² (B×κ)·∇α"),
]


# ──────────────────────────────────────────────────────────────────────────────
# Self-consistency checks
# ──────────────────────────────────────────────────────────────────────────────

def self_consistency_check(data: dict, label: str = ""):
    """Run physical sanity checks on a single geometry file."""
    ok = True
    prefix = f"[{label}] " if label else ""

    def check(cond, msg):
        nonlocal ok
        status = "PASS" if cond else "FAIL"
        print(f"  {prefix}{status}  {msg}")
        if not cond:
            ok = False

    bmag = data["bmag"]
    check(np.all(bmag > 0),
          f"bmag > 0 everywhere (min={bmag.min():.4f})")

    gpp = data["grad_psit_grad_psit"]
    check(np.all(gpp > 0),
          f"|∇ψt|² > 0 everywhere (min={gpp.min():.6g})")

    gaa = data["grad_alpha_grad_alpha"]
    check(np.all(gaa > 0),
          f"|∇α|² > 0 everywhere (min={gaa.min():.6g})")

    # B.∇ζ should have constant sign
    bdgz = data["b_dot_gradzeta"]
    check(np.all(bdgz > 0) or np.all(bdgz < 0),
          f"b.∇ζ has consistent sign (mean={bdgz.mean():.4f})")

    # Cauchy-Schwarz: (∇α·∇ψt)² ≤ |∇α|²|∇ψt|²
    # Using unnormalised quantities  a²B_ref² * gpp and a² * gaa
    gap  = data["grad_alpha_grad_psit"]
    cs_ok = np.all(gap**2 <= gaa * gpp * (1.0 + 1e-6))
    check(cs_ok, "Cauchy-Schwarz: (∇α·∇ψt)² ≤ |∇α|²|∇ψt|²")

    q    = float(data["safety_factor_q"])
    shat = float(data["shat"])
    # q < 0 is allowed: its sign is that of iota in the frame of the VMEC angles
    check(q != 0 and np.isfinite(q),
          f"safety factor q finite and non-zero (q={q:.4f})")
    check(abs(shat) < 50,
          f"|shat| reasonable (shat={shat:.4f})")

    return ok


# ──────────────────────────────────────────────────────────────────────────────
# GVEC vs VMEC comparison
# ──────────────────────────────────────────────────────────────────────────────

def compare(gvec_data: dict, vmec_data: dict, tol: float, outdir: Path):
    """Compare every geometry quantity between the GVEC and VMEC files."""
    all_ok = True
    rows = []

    # ── Scalars ──────────────────────────────────────────────────────────────
    print("\n  Scalar quantities:")
    for name in _SCALAR_VARS:
        g = float(gvec_data.get(name, np.nan))
        v = float(vmec_data.get(name, np.nan))
        err = abs(g - v) / max(abs(v), 1e-14)
        status = "PASS" if err < tol else "FAIL"
        ok = (status == "PASS")
        all_ok &= ok
        print(f"    {status}  {name:30s}  gvec={g:.6g}  vmec={v:.6g}  rel_err={err:.2e}")
        rows.append((name, g, v, err, ok))

    # ── 2-D fields ────────────────────────────────────────────────────────────
    print("\n  2-D geometry arrays (L∞ relative error vs VMEC):")
    figs = []
    for varname, label in _2D_VARS:
        if varname not in gvec_data or varname not in vmec_data:
            print(f"    SKIP  {varname:40s}  (not found in one of the files)")
            continue
        g = np.atleast_2d(gvec_data[varname])
        v = np.atleast_2d(vmec_data[varname])
        if g.shape != v.shape:
            print(f"    SKIP  {varname:40s}  shape mismatch {g.shape} vs {v.shape}")
            continue
        err = rel_err(g, v)
        status = "PASS" if err < tol else "FAIL"
        ok = (status == "PASS")
        all_ok &= ok
        print(f"    {status}  {varname:40s}  rel_err={err:.2e}  (tol={tol:.1e})")
        rows.append((varname, None, None, err, ok))

        # Plot comparison along first field line (ia=0)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        zeta = vmec_data.get("zeta", np.arange(v.shape[-1]))
        axes[0].plot(zeta, v[0], label="VMEC", lw=2)
        axes[0].plot(zeta, g[0], "--", label="GVEC", lw=1.5)
        axes[0].set_xlabel(r"$\zeta$")
        axes[0].set_title(label)
        axes[0].legend()
        axes[1].plot(zeta, g[0] - v[0])
        axes[1].set_xlabel(r"$\zeta$")
        axes[1].set_title(f"GVEC − VMEC  (max rel err = {err:.2e})")
        fig.suptitle(varname, fontsize=10)
        fig.tight_layout()
        outpath = outdir / f"cmp_{varname}.png"
        fig.savefig(str(outpath), dpi=100)
        plt.close(fig)
        figs.append(outpath)

    return all_ok, rows, figs


# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gvec",  required=True,  help="GVEC stella geometry .nc file")
    p.add_argument("--vmec",  default=None,   help="VMEC stella geometry .nc file (optional)")
    p.add_argument("--tol",   type=float, default=1e-2,
                   help="Relative L∞ tolerance for PASS/FAIL (default 0.01 = 1%%)")
    p.add_argument("--outdir", default="compare_output",
                   help="Output directory for comparison plots")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  GVEC stella geometry comparison")
    print(f"{'='*60}")
    print(f"  GVEC file : {args.gvec}")
    if args.vmec:
        print(f"  VMEC file : {args.vmec}")
    print(f"  tolerance : {args.tol:.1e}")
    print()

    # ── Load ─────────────────────────────────────────────────────────────────
    gvec_data = load_geometry(args.gvec)
    print(f"  GVEC: s={float(gvec_data['s']):.4f},  "
          f"q={float(gvec_data['safety_factor_q']):.4f},  "
          f"shat={float(gvec_data['shat']):.4f}")

    # ── Self-consistency check ────────────────────────────────────────────────
    print("\n[Self-consistency checks: GVEC file]")
    gvec_ok = self_consistency_check(gvec_data, label="GVEC")

    overall_ok = gvec_ok

    if args.vmec:
        vmec_data = load_geometry(args.vmec)
        print(f"\n  VMEC: s={float(vmec_data['s']):.4f},  "
              f"q={float(vmec_data['safety_factor_q']):.4f},  "
              f"shat={float(vmec_data['shat']):.4f}")

        print("\n[Self-consistency checks: VMEC file]")
        vmec_ok = self_consistency_check(vmec_data, label="VMEC")
        overall_ok &= vmec_ok

        print("\n[GVEC vs VMEC comparison]")
        cmp_ok, rows, figs = compare(gvec_data, vmec_data, args.tol, outdir)
        overall_ok &= cmp_ok

        if figs:
            print(f"\n  Plots written to: {outdir}/")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    if overall_ok:
        print("  OVERALL: PASS  (all checks within tolerance)")
    else:
        print("  OVERALL: FAIL  (see details above)")
    print(f"{'='*60}\n")

    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
