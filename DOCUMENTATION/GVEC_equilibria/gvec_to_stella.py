#!/usr/bin/env python3
# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
r"""gvec_to_stella.py
====================
Post-processing script: convert a GVEC equilibrium state file into a geometry
netCDF file that stella can read directly (via ``geometry_option = 'gvec'``).

The script evaluates all geometry arrays that stella needs – |B|, gradient
dot-products, drift terms, etc. – on the (alpha, zeta) field-line grid stella
will use, and writes them to ``<output>.nc``.  The Fortran reader
``geometry_gvec_read_netCDF_file.f90`` then reads this file at run-time;
no linking against GVEC is required.

Usage
-----
::

    python gvec_to_stella.py  STATE_FILE  [options]

Key options
~~~~~~~~~~~
``--s``                 Normalised toroidal flux of the desired surface (0–1).
``--nalpha``            Number of field lines.
``--alpha0``            First field-line label (radians).
``--nzgrid``            Half number of grid points along z  (total: 2*nzgrid+1).
``--nzgrid_file``       Half number of points written to file (>= nzgrid; used
                        as a refinement factor for arc-length interpolation in
                        stella).  Defaults to nzgrid.
``--zeta_center``       Centre of the toroidal domain (radians).
``--nfield_periods``    Number of field periods to include (<= nfp from file).
                        Use -1 (default) to use all nfp field periods.
``--output``            Output file name (default: gvec_stella_geometry.nc).

Notes
-----
*   All output quantities follow stella's VMEC normalisations:
    - L_reference = a_minor  (metres)
    - B_reference = 2 |phi_edge/(2pi)| / a_minor^2  (Tesla)
    - grad_alpha_grad_alpha = a^2 |grad alpha|^2
    - etc.  (see geometry_vmec_read_netCDF_file.f90 comments for definitions)

*   For arc-length coordinates in stella (``zed_equal_arc = .true.``) run
    this script with ``--nzgrid_file N*M`` where M is the refinement factor
    and N is stella's nzgrid.  Stella's Fortran code will then do the
    arc-length interpolation at run-time.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import netCDF4 as nc

try:
    import gvec
    from gvec.core.compute import EvaluationsPEST, compute
except ImportError:
    sys.exit(
        "ERROR: pygvec is not installed.  "
        "Install it from the gvec repository: pip install -e gvec/"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Helper: dot product over the xyz (axis-0) dimension
# ──────────────────────────────────────────────────────────────────────────────

def _dot(A, B):
    """Dot product over axis 0 (xyz); A, B have shape (3, ...)."""
    return np.einsum("i...,i...->...", A, B)


def _cross(A, B):
    """Cross product over axis 0 (xyz); A, B have shape (3, ...) → (3, ...)."""
    return np.cross(A, B, axis=0)


def _get_array(ev, name):
    """Extract a quantity from the xarray Dataset, squeezing the rad dimension."""
    da = ev[name]
    arr = da.values
    # Scalar fields may have dims (rad, pol, tor) → squeeze rad
    if "rad" in da.dims:
        rad_ax = da.dims.index("rad")
        arr = np.squeeze(arr, axis=rad_ax)
    return arr  # shape (nalpha, nz) or (nalpha, nz, 3) etc.


def _get_vec(ev, name):
    """Extract a Cartesian vector (3, nalpha, nz) from the Dataset."""
    da = ev[name]
    # Move xyz to first axis
    arr = da.transpose("xyz", ...).values
    # Drop rad axis (size 1) from positions 1+ if present
    while arr.ndim > 3:
        arr = arr[:, 0, :, :]  # remove leading rad dim
    return arr  # (3, nalpha, nz)


# ──────────────────────────────────────────────────────────────────────────────
# Main converter
# ──────────────────────────────────────────────────────────────────────────────

def _resolve_parameter_file(state_file, parameter_file=None):
    """Locate the GVEC parameter .ini that goes with a state .dat file.

    ``gvec.load_state`` needs BOTH the parameter file and the state file.
    If the caller did not supply one, look next to the state file for a
    ``parameter*.ini``, preferring a name containing "final".
    """
    state_file = Path(state_file)
    if parameter_file is not None:
        parameter_file = Path(parameter_file)
        if not parameter_file.is_file():
            sys.exit(f"ERROR: parameter file not found: {parameter_file}")
        return parameter_file

    candidates = sorted(state_file.parent.glob("parameter*.ini"))
    if not candidates:
        candidates = sorted(state_file.parent.glob("*.ini"))
    if not candidates:
        sys.exit(
            f"ERROR: no GVEC parameter (.ini) file found next to {state_file}.\n"
            "       Pass one explicitly with --parameter-file."
        )
    for c in candidates:
        if "final" in c.name:
            return c
    return candidates[0]


def gvec_to_stella(
    state_file: str | Path,
    parameter_file: str | Path | None = None,
    a_minor_override: float | None = None,
    flip: str = "auto",
    s_input: float = 0.5,
    nalpha: int = 1,
    alpha0: float = 0.0,
    nzgrid: int = 128,
    nzgrid_file: int = -1,          # -1 → same as nzgrid
    zeta_center: float = 0.0,
    nfield_periods: int = -1,        # -1 → use all nfp
    output_file: str | Path | None = None,
    verbose: bool = True,
) -> Path:
    """Convert a GVEC state file to a stella geometry netCDF file.

    Parameters
    ----------
    state_file : path-like
        Path to the GVEC state (restart) file, e.g. ``PROJ_State_00001.dat``.
    s_input : float
        Normalised toroidal flux (0–1) of the desired flux surface.
    nalpha : int
        Number of field lines.
    alpha0 : float
        Label of the first field line (radians).
    nzgrid : int
        Half number of z-grid points that *stella* will use.
    nzgrid_file : int
        Half number of z-grid points written to the file.  Must satisfy
        ``nzgrid_file >= nzgrid``.  Use ``nzgrid_file = nzgrid * M`` for
        arc-length refinement factor M.  Defaults to ``nzgrid``.
    zeta_center : float
        Centre of the toroidal domain (radians).
    nfield_periods : int
        Number of field periods to simulate (≤ nfp from state).
        -1 means use all field periods (nfp).
    output_file : path-like or None
        Output netCDF filename.  Defaults to ``gvec_stella_geometry.nc``.
    verbose : bool
        Print progress messages.

    Returns
    -------
    Path
        Path to the written netCDF file.
    """
    state_file = Path(state_file)
    if nzgrid_file < 0:
        nzgrid_file = nzgrid
    if nzgrid_file < nzgrid:
        raise ValueError(f"nzgrid_file={nzgrid_file} must be >= nzgrid={nzgrid}")
    if output_file is None:
        output_file = Path("gvec_stella_geometry.nc")
    output_file = Path(output_file)

    mu0 = 4.0 * np.pi * 1.0e-7   # SI

    # ── 1. Load state ─────────────────────────────────────────────────────────
    if verbose:
        print(f"Loading GVEC state: {state_file}")
    parameter_file = _resolve_parameter_file(state_file, parameter_file)
    if verbose:
        print(f"  using GVEC parameter file: {parameter_file}")
    state = gvec.load_state(str(parameter_file), str(state_file))
    nfp = state.nfp
    nfp_sim = nfield_periods if nfield_periods > 0 else nfp

    # ── 2. Radial quantities ──────────────────────────────────────────────────
    # GVEC uses rho = sqrt(Phi/Phi_edge) ∈ [0,1] as radial coordinate,
    # where Phi is the toroidal flux.
    # Stella's VMEC uses  s = Phi/Phi_edge = rho² ∈ [0,1].
    rho0 = float(np.sqrt(s_input))
    if verbose:
        print(f"  surface: s = {s_input:.4f}  (rho = {rho0:.4f})")

    iota0       = float(state.evaluate_profile("iota", np.array([rho0]))[0])
    diota_drho  = float(state.evaluate_profile("iota", np.array([rho0]), deriv=1)[0])
    # NOTE: GVEC's "Phi" profile is the toroidal flux PER RADIAN, i.e. it is
    # already Phi_edge/(2*pi) - not the total flux. (Cross-check: the .ini for
    # this case sets phiedge = -65.7029..., while evaluate_profile returns
    # -10.4569... = -65.7029.../(2*pi).) GVEC's own gist.py exporter relies on
    # the same convention when it forms Bref = 2*Phi_edge/a**2.
    phi_edge_over_2pi = float(state.evaluate_profile("Phi", np.array([1.0]))[0])
    Phi_edge    = phi_edge_over_2pi * 2.0 * np.pi                            # Wb (total)
    dp_dr       = float(state.evaluate_profile("p",    np.array([rho0]), deriv=1)[0])

    # ── Toroidal/poloidal handedness ──────────────────────────────────────────
    # The formulas of stella's VMEC geometry hold in the frame of the VMEC
    # angles (theta, zeta), whose Jacobian is negative (signgs = -1), with the
    # toroidal flux of the right-handed frame, psi_t = signgs * phi_VMEC, as the
    # radial coordinate (geometry_vmec_read_netCDF_file.f90). GVEC's logical
    # coordinates (rho, theta, zeta) are right-handed (positive Jacobian), so one
    # angle is reversed to obtain a frame of the same handedness as VMEC's; with
    # the cylindrical map of GVEC (hmap = 1), reversing zeta gives the VMEC angles
    # of the same equilibrium (e.g. after pygvec convert-wout). Reversing an angle
    # changes the sign of iota and of the Jacobian, but not the toroidal flux
    # Phi' = B^zeta * Jac, so psi_t is the toroidal flux of GVEC in both frames.
    #
    # A frame with a positive Jacobian (no flip, or both angles flipped) gives the
    # same metric coefficients and drifts, but the wrong sign of the Clebsch
    # factor B = sign_torflux grad(psi_t) x grad(alpha) for stella, which reverses
    # the drive relative to the drifts, so it is not offered.
    _ev_j = state.evaluate("Jac", rho=rho0, theta=np.array([0.0]), zeta=np.array([0.0])).squeeze()
    if float(_ev_j.Jac) <= 0:
        sys.exit("ERROR: the Jacobian of the GVEC state is not positive.")
    if flip == "auto":
        flip = "tor"
    if flip not in ("tor", "pol"):
        sys.exit(f"ERROR: flip must be tor/pol/auto, got '{flip}'")

    flip_tor = flip == "tor"
    flip_pol = flip == "pol"

    # Reversing one angle reverses the sign of iota.
    iota0 = -iota0
    diota_drho = -diota_drho
    if verbose:
        print(f"  flip = '{flip}'")

    q0   = 1.0 / iota0
    # shat = (rho/q)(dq/drho) = -rho * (diota/drho) / iota
    shat = -rho0 * diota_drho / iota0
    # d(iota)/d(s_vmec) where s_vmec = rho²
    diota_ds = diota_drho / (2.0 * rho0)
    # dp/d(s_vmec)
    dp_ds = dp_dr / (2.0 * rho0)

    # Reference quantities (same as stella's VMEC convention)
    # Effective minor radius (GVEC's "r_minor"); this is a global constant of
    # the equilibrium, so it can be evaluated on any grid. Matches the
    # convention used by GVEC's own gist.py exporter (Lref = r_minor).
    _ev_a = gvec.evaluate(state, rho=1.0, theta="int", zeta="int")
    state.compute(_ev_a, "r_minor")
    a_minor = float(np.asarray(_ev_a["r_minor"]).ravel()[0])

    # GVEC's r_minor is a VOLUME-based effective radius, sqrt(V/(pi*L_axis)),
    # using the true length of the magnetic axis. VMEC's Aminor_p is
    # AREA-based and effectively uses 2*pi*Rmajor_p instead of L_axis. For a
    # non-planar axis (W7-X: L_axis = 142.0 m vs 2*pi*R0 = 134.1 m) the two
    # differ by ~3%, which propagates to ~6% in B_ref (B_ref ~ 1/a^2) and so
    # rescales every normalised geometry array.
    #
    # This is a definitional difference, not an error in either code - but a
    # like-for-like benchmark against stella's VMEC branch needs both runs on
    # the SAME normalisation. Pass --a-minor <VMEC Aminor_p> to do that.
    if a_minor_override is not None:
        if verbose:
            print(f"  overriding a_minor: {a_minor:.6f} m (GVEC r_minor) "
                  f"-> {a_minor_override:.6f} m (requested)")
        a_minor = float(a_minor_override)
    B_ref   = 2.0 * abs(phi_edge_over_2pi) / a_minor**2   # Tesla
    L_ref   = a_minor                                       # metres
    sign_torflux = int(np.sign(Phi_edge))

    if verbose:
        print(f"  q = {q0:.4f},  shat = {shat:.4f}")
        print(f"  L_ref = {L_ref:.4f} m,  B_ref = {B_ref:.4f} T")
        print(f"  sign(torflux) = {sign_torflux},  nfp = {nfp},  nfp_sim = {nfp_sim}")

    # ── 3. Grid ───────────────────────────────────────────────────────────────
    nz_file = 2 * nzgrid_file + 1
    alpha   = alpha0 + np.arange(nalpha) * 2.0 * np.pi / nalpha   # (nalpha,)
    # zeta grid: same formula as stella's VMEC reader
    zeta = (zeta_center
            + np.pi * np.arange(-nzgrid_file, nzgrid_file + 1)
            * nfp_sim / (nfp * nzgrid_file))                        # (nz_file,)

    # theta_P on each field line: theta_P(alpha, zeta) = alpha + iota * zeta
    theta_P_2d = alpha[:, np.newaxis] + iota0 * zeta[np.newaxis, :]  # (nalpha, nz_file)

    # <zeta>/<theta_P_2d> above are in the FLIPPED (stella-facing) frame. GVEC
    # must be evaluated in its own frame, so undo the flip on the sample points.
    zeta_gvec    = -zeta       if flip_tor else zeta
    theta_P_gvec = -theta_P_2d if flip_pol else theta_P_2d

    # ── 4. Evaluate geometry at PEST grid points ───────────────────────────────
    if verbose:
        print(f"  Setting up PEST grid ({nalpha} x {nz_file} points)...")

    # EvaluationsPEST performs the Newton solve theta(theta_P, zeta) internally.
    ev = EvaluationsPEST(
        rho      = rho0,          # scalar → rad dim of size 1
        theta_P  = theta_P_gvec,  # 2-D array (nalpha, nz_file) → dims ("pol", "tor")
        zeta     = zeta_gvec,     # 1-D array (nz_file,)       → dim  ("tor",)
        state    = state,
    )

    if verbose:
        print("  Computing geometry quantities...")

    # NOTE: compute() adds the quantities to <ev> IN PLACE and returns None,
    # so its result must not be assigned back to <ev>.
    compute(
        ev,
        "pos",   # position vector; "xyz" only registers the (x,y,z) coord labels
        "grad_rho",
        "grad_theta_P",
        "grad_zeta",
        "B",
        "B_contra_t",
        "B_contra_z",
        "mod_B",
        "grad_mod_B",
        "g_tt",
        "g_tz",
        "g_zz",
        "Jac",
        "dLA_dr",
        "theta",
        state=state,
    )

    # ── 5. Extract arrays ─────────────────────────────────────────────────────
    xyz_arr      = _get_vec(ev, "pos")        # (3, nalpha, nz_file)  Cartesian pos
    grad_rho     = _get_vec(ev, "grad_rho")
    grad_theta_P = _get_vec(ev, "grad_theta_P")
    grad_zeta_v  = _get_vec(ev, "grad_zeta")
    B_cart       = _get_vec(ev, "B")
    grad_modB    = _get_vec(ev, "grad_mod_B")

    modB        = _get_array(ev, "mod_B")     # (nalpha, nz_file)
    B_ct        = _get_array(ev, "B_contra_t")
    B_cz        = _get_array(ev, "B_contra_z")
    g_tt        = _get_array(ev, "g_tt")
    g_tz        = _get_array(ev, "g_tz")
    g_zz        = _get_array(ev, "g_zz")
    dLA_dr_arr  = _get_array(ev, "dLA_dr")
    theta_arr   = _get_array(ev, "theta")     # logical poloidal angle

    # ── Express the evaluated fields in the FLIPPED frame ─────────────────────
    # Quantities carrying a free theta or zeta index change sign; true vectors
    # and scalars (pos, grad_rho, B, grad|B|, |B|) and the angle-angle metric
    # components (both indices flip, so the sign cancels) do not.
    if flip_tor:
        grad_zeta_v = -grad_zeta_v      # grad zeta_new = -grad zeta_gvec
        B_cz        = -B_cz             # B^zeta
    if flip_pol:
        grad_theta_P = -grad_theta_P    # grad theta_new = -grad theta_gvec
        B_ct         = -B_ct            # B^theta
        theta_arr    = -theta_arr       # the logical poloidal angle itself
        dLA_dr_arr   = -dLA_dr_arr      # lambda flips with theta

    # Cylindrical radius R = sqrt(x²+y²)
    R = np.sqrt(xyz_arr[0]**2 + xyz_arr[1]**2)   # (nalpha, nz_file)

    # ── 6. Derive stella geometry arrays ─────────────────────────────────────
    # ∇s  = 2ρ ∇ρ          (s_vmec = ρ²)
    # ∇ψt = (Phi_edge/2π) ∇s = 2ρ (Phi_edge/2π) ∇ρ
    grad_s    = 2.0 * rho0 * grad_rho
    grad_psit = phi_edge_over_2pi * grad_s         # (3, nalpha, nz_file)

    # ∇α = ∇θ_P - ι ∇ζ - (ζ − ζ_c) (dι/dρ) ∇ρ
    # (centred at zeta_center, consistent with stella's VMEC formula)
    # zeta - zeta_center has shape (1, nz_file) → broadcast to (nalpha, nz_file)
    zeta_minus_c = (zeta - zeta_center)[np.newaxis, :]   # (1, nz_file)
    grad_alpha = (grad_theta_P
                  - iota0       * grad_zeta_v
                  - zeta_minus_c * diota_drho * grad_rho)  # (3, nalpha, nz_file)

    # ── Scalar normalised quantities ──────────────────────────────────────────
    bmag = modB / B_ref                            # B̃ = |B|/Bref

    # b · ∇ζ = (a/|B|) B·∇ζ = L_ref * B^ζ / |B|
    b_dot_gradzeta = L_ref * B_cz / modB

    # Grad dot products (stella normalisation – see calculate_vmec_geometry comments)
    gag = L_ref**2 * _dot(grad_alpha, grad_alpha)         # a² |∇α|²
    gap = (1.0 / B_ref) * _dot(grad_alpha, grad_psit)     # ∇α·∇ψt / Bref
    gpp = 1.0 / (L_ref**2 * B_ref**2) * _dot(grad_psit, grad_psit)  # |∇ψt|²/(a²Bref²)

    # Intermediate: ∇ζ·∇α and ∇ζ·∇ψt  (unnormalised analogue for gds23-26)
    gza = L_ref       * _dot(grad_zeta_v, grad_alpha)
    gzp = 1.0 / (L_ref * B_ref) * _dot(grad_zeta_v, grad_psit)
    gta = L_ref       * _dot(grad_theta_P, grad_alpha)
    gtp = 1.0 / (L_ref * B_ref) * _dot(grad_theta_P, grad_psit)

    bmag2 = bmag**2
    gds23 = -1.0 / bmag2 * (gza * gap - gzp * gag)
    gds24 = -1.0 / bmag2 * (gza * gpp - gzp * gap)
    gds25 = -sign_torflux / bmag2 * (gta * gap - gtp * gag)
    gds26 = -0.5 * sign_torflux / bmag2 * (gta * gpp - gtp * gap)

    # Drift terms: (B×∇B)·∇α  and  (B×∇B)·∇ψt
    BcGB     = _cross(B_cart, grad_modB)
    BcGBdGA  = _dot(BcGB, grad_alpha)
    BcGBdGP  = _dot(BcGB, grad_psit)

    # Drift terms: (B×∇s)·∇α  for curvature including pressure gradient
    BcGs     = _cross(B_cart, grad_s)
    BcGsdGA  = _dot(BcGs, grad_alpha)
    BcGsdGP  = _dot(BcGs, grad_psit)

    # B_times_gradB_dot_gradalpha = a²Bref/|B|³ (B×∇B)·∇α
    Btg_dGA  = L_ref**2 * B_ref / modB**3 * BcGBdGA
    Btg_dGP  = 1.0 / modB**3 * BcGBdGP

    # B_times_kappa_dot_gradalpha includes the pressure gradient term
    # B_times_kappa = B_times_gradB + (a²Bref μ₀/B⁴)(dp/ds)(B×∇s)
    Btk_dGA = Btg_dGA + L_ref**2 * B_ref * mu0 / modB**4 * dp_ds * BcGsdGA
    Btk_dGP = Btg_dGP + mu0 / modB**4 * dp_ds * BcGsdGP

    # Covariant B components: B_θ = g_θθ B^θ + g_θζ B^ζ, B_ζ = g_θζ B^θ + g_ζζ B^ζ
    B_sub_theta = g_tt * B_ct + g_tz * B_cz
    B_sub_zeta  = g_tz * B_ct + g_zz * B_cz

    # R²/B² terms (for momentum flux in stella)
    R2oB2              = R**2 / bmag2
    gza_R2oB2          = gza * R2oB2
    gzp_R2oB2          = gzp * R2oB2
    b_dot_gradzeta_RR  = b_dot_gradzeta * R**2

    # psit_displacement_fac = (∂λ/∂s_vmec − (ζ−ζ_c) dι/ds_vmec) ψt_edge
    # = (dLA_dr/(2ρ) − (ζ−ζ_c) diota_drho/(2ρ)) phi_edge_over_2pi
    psit_disp = (dLA_dr_arr - zeta_minus_c * diota_drho) * phi_edge_over_2pi / (2.0 * rho0)

    # grho = a |∇ρ|  (= sqrt(grad_psit_grad_psit) / rhotor, confirmed by VMEC code)
    grho = L_ref * np.sqrt(_dot(grad_rho, grad_rho))

    # ── 7. Compute dzdzeta (conversion factor z = zeta*P) ─────────────────────
    # P = nfp / nfp_sim  so that z is compressed to [-π, π]
    dzdzeta = float(nfp) / nfp_sim

    # ── 8. Write netCDF ───────────────────────────────────────────────────────
    if verbose:
        print(f"  Writing: {output_file}")

    if output_file.exists():
        output_file.unlink()

    ds = nc.Dataset(str(output_file), "w")
    try:
        # ── dimensions ──
        ds.createDimension("nalpha", nalpha)
        ds.createDimension("nz_file", nz_file)   # 2*nzgrid_file + 1

        # ── integer scalars ──
        for name, val in [
            ("nfp",          nfp),
            ("sign_torflux", sign_torflux),
            ("nalpha",       nalpha),
            ("nzgrid",       nzgrid),       # stella's nzgrid (for validation)
            ("nzgrid_file",  nzgrid_file),  # half-points in this file
            ("nfp_sim",      nfp_sim),
        ]:
            v = ds.createVariable(name, "i4")
            v.assignValue(int(val))

        # ── real scalars ──
        for name, val in [
            ("s",                s_input),
            ("safety_factor_q",  q0),
            ("shat",             shat),
            ("L_reference",      L_ref),
            ("B_reference",      B_ref),
            ("alpha0",           alpha0),
            ("zeta_center",      zeta_center),
            ("dzdzeta",          dzdzeta),
        ]:
            v = ds.createVariable(name, "f8")
            v.assignValue(float(val))

        # ── 1-D arrays ──
        v = ds.createVariable("alpha", "f8", ("nalpha",))
        v[:] = alpha
        v.long_name = "field line label alpha = theta_P - iota*zeta"

        v = ds.createVariable("zeta", "f8", ("nz_file",))
        v[:] = zeta
        v.long_name = "toroidal angle grid"

        # ── 2-D arrays (nalpha, nz_file) ──
        _2d_vars = {
            "theta":                      (theta_arr,        "logical poloidal angle"),
            "grho":                       (grho,             "a |grad rho|"),
            "bmag":                       (bmag,             "B/Bref"),
            "b_dot_gradzeta":             (b_dot_gradzeta,   "(a/B) B.grad(zeta)"),
            "B_sub_theta":                (B_sub_theta,      "covariant poloidal B"),
            "B_sub_zeta":                 (B_sub_zeta,       "covariant toroidal B"),
            "grad_alpha_grad_alpha":      (gag,              "a^2 |grad alpha|^2"),
            "grad_alpha_grad_psit":       (gap,              "(1/Bref) grad_alpha.grad_psit"),
            "grad_psit_grad_psit":        (gpp,              "|grad psit|^2 / (a^2 Bref^2)"),
            "gds23":                      (gds23,            "gds23_psitalpha"),
            "gds24":                      (gds24,            "gds24_psitalpha"),
            "gds25":                      (gds25,            "gds25_psitalpha"),
            "gds26":                      (gds26,            "gds26_psitalpha"),
            "gbdrift_alpha":              (Btg_dGA,          "a^2 Bref/B^3 (BxgradB).grad_alpha"),
            "B_times_gradB_dot_gradpsit": (Btg_dGP,         "1/B^3 (BxgradB).grad_psit"),
            "gbdrift_kappa_alpha":        (Btk_dGA,         "a^2 Bref/B^2 (Bxkappa).grad_alpha"),
            "B_times_kappa_dot_gradpsit": (Btk_dGP,         "1/B^2 (Bxkappa).grad_psit"),
            "psit_displacement_fac":      (psit_disp,        "psit displacement factor"),
            "gradzeta_gradalpha_R2overB2":(gza_R2oB2,        "(a*R^2/B^2) grad_zeta.grad_alpha"),
            "gradzeta_gradpsit_R2overB2": (gzp_R2oB2,        "(R^2/aB^2) grad_zeta.grad_psit"),
            "b_dot_gradzeta_RR":          (b_dot_gradzeta_RR,"b.grad_zeta * R^2"),
        }
        for name, (arr, longname) in _2d_vars.items():
            v = ds.createVariable(name, "f8", ("nalpha", "nz_file"))
            v[:, :] = arr
            v.long_name = longname

        # ── metadata ──
        ds.state_file    = str(state_file)
        ds.s_input       = s_input
        ds.nfp           = nfp
        ds.description   = (
            "GVEC equilibrium geometry for stella gyrokinetic code. "
            "Generated by gvec_to_stella.py."
        )

    finally:
        ds.close()

    if verbose:
        print(f"Done. Output: {output_file}")

    return output_file


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("state_file", help="Path to GVEC state file")
    p.add_argument("--flip", default="auto",
                   choices=["auto", "tor", "pol"],
                   help="reverse the toroidal (tor) or poloidal (pol) angle, so that the "
                        "frame has the handedness of the VMEC angles that stella's VMEC "
                        "conventions assume (default: auto = tor, which gives the VMEC "
                        "angles for the cylindrical map of GVEC)")
    p.add_argument("--a-minor", "--a_minor", type=float, default=None,
                   help="override the reference minor radius in metres (e.g. VMEC's "
                        "Aminor_p) so a GVEC export shares a VMEC run's normalisation")
    p.add_argument("--parameter-file", "--parameter_file", default=None,
                   help="GVEC parameter .ini file (default: auto-detect next to the state file)")
    p.add_argument("--s",              type=float, default=0.5,
                   help="Normalised toroidal flux of the surface (default: 0.5)")
    p.add_argument("--nalpha",         type=int,   default=1,
                   help="Number of field lines (default: 1)")
    p.add_argument("--alpha0",         type=float, default=0.0,
                   help="Label of first field line in radians (default: 0)")
    p.add_argument("--nzgrid",         type=int,   default=128,
                   help="Half number of z-grid points for stella (default: 128)")
    p.add_argument("--nzgrid_file",    type=int,   default=-1,
                   help="Half number of z-grid points in file (default: same as nzgrid). "
                        "Set to nzgrid*M to enable arc-length refinement factor M in stella.")
    p.add_argument("--zeta_center",    type=float, default=0.0,
                   help="Centre of the toroidal domain in radians (default: 0)")
    p.add_argument("--nfield_periods", type=int,   default=-1,
                   help="Number of field periods to include (-1 = all, default: -1)")
    p.add_argument("--output", "-o",   default=None,
                   help="Output netCDF filename (default: gvec_stella_geometry.nc)")
    p.add_argument("--quiet",  "-q",   action="store_true",
                   help="Suppress progress output")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    gvec_to_stella(
        state_file     = args.state_file,
        parameter_file = args.parameter_file,
        a_minor_override = args.a_minor,
        flip           = args.flip,
        s_input        = args.s,
        nalpha         = args.nalpha,
        alpha0         = args.alpha0,
        nzgrid         = args.nzgrid,
        nzgrid_file    = args.nzgrid_file,
        zeta_center    = args.zeta_center,
        nfield_periods = args.nfield_periods,
        output_file    = args.output,
        verbose        = not args.quiet,
    )
