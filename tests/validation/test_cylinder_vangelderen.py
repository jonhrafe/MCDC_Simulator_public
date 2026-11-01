#!/usr/bin/env python3
"""
Analytical validation: restricted diffusion in a single infinite cylinder.

This is the CYL-chain template. It runs tests/accuracy/cylinder_restricted.conf
(one native-list infinite z-cylinder, radius 5 um, walkers seeded intra,
IMPERMEABLE) and validates the restricted-diffusion signal two ways:

(a) ANISOTROPY SANITY (always on, no dmipy needed).
    A z-oriented cylinder is anisotropic: a gradient ALONG z probes free axial
    diffusion (S/S0 ~ exp(-bD), strong attenuation) while a PERPENDICULAR
    gradient probes restricted diffusion (weak attenuation, high S/S0). So at
    every shell the perpendicular signal must sit clearly ABOVE both the parallel
    signal and the free curve. This is a cheap, robust physics gate.

(b) vs the frozen van Gelderen curve (skips if the curve is absent).
    Compares the full angular signal pattern per shell to the dmipy
    C4CylinderGaussianPhaseApproximation curve frozen at tests/accuracy/analytic/
    (see tests/tools/gen_analytic.py). Per-shell NRMSE within a statistical band.
    Valid only at small b (the GPA is an infinite impermeable cylinder).

Tolerances are statistical: fixed-step MC has noise ~1/sqrt(N) plus an O(1/T)
finite-step bias that grows with b -- this validates the physics, not
bit-exactness. Runs as pytest (CPU + GPU-self-skip) or standalone:

  python3 tests/validation/test_cylinder_vangelderen.py --gpu 0 --n 20000
"""
import argparse
import math
import os
import shutil
import sys
import tempfile
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcdc_testkit import io, metrics, oracles, paths, run, scheme  # noqa: E402

import numpy as np  # noqa: E402

CONF = paths.accuracy("cylinder_restricted.conf")
ANALYTIC_NAME = "cylinder_r5um_d0p6"
RADIUS = 5e-6         # single_cylinder.list radius (matches ANALYTIC_NAME r5um)
PAR_MAX_DEG = 30.0    # angle-to-z <= 30 deg -> gradient ~along axis (free axial diffusion)
PERP_MIN_DEG = 60.0   # angle-to-z >= 60 deg -> gradient ~perpendicular (restricted)

Result = namedtuple("Result", "ss0 sc D n skipped")


def run_cylinder(binary, n=20000, t=None, gpu=0, timeout=600):
    """Run the single-cylinder case; return a :class:`Result` (normalised S/S0)."""
    conf_text = open(CONF).read()
    N = int(n) if n else int(io.read_conf_value(conf_text, "N"))
    D = float(io.read_conf_value(conf_text, "diffusivity"))
    sc = scheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))
    overrides = {"N": N}
    if t:
        overrides["T"] = int(t)
    work = tempfile.mkdtemp(prefix="mcdc_cyl_")
    try:
        prefix = run.run_case(binary, conf_text, work, "cyl", paths.REPO_ROOT,
                              overrides=overrides, gpu=gpu, timeout=timeout)
        if prefix is None:
            return Result(None, None, None, 0, True)
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=N)
        if dwi.size != sc.n:
            raise AssertionError(f"DWI entries ({dwi.size}) != scheme rows ({sc.n})")
        return Result(metrics.normalize_by_b0(dwi, sc.b0_mask), sc, D, sc.n, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def anisotropy_rows(r):
    """Per-shell (b>0) mean S/S0 for parallel-ish vs perpendicular-ish gradients + free exp(-bD)."""
    sc = r.sc
    nb = ~sc.b0_mask
    ang = scheme.angle_to_axis_deg(sc.gdirs)
    b_mm = np.round(sc.bvals / 1e6, 0)
    rows = []
    for b in sorted(set(b_mm[nb])):
        if b <= 0:
            continue
        m = nb & (b_mm == b)
        par = r.ss0[m & (ang <= PAR_MAX_DEG)]
        per = r.ss0[m & (ang >= PERP_MIN_DEG)]
        rows.append(dict(b=float(b), n_par=int(par.size), n_per=int(per.size),
                         par=float(par.mean()) if par.size else float("nan"),
                         per=float(per.mean()) if per.size else float("nan"),
                         free=math.exp(-(b * 1e6) * r.D)))
    return rows


def check_anisotropy(rows):
    """(ok, messages). Require >=2 shells (both dir-groups populated) where the
    perpendicular restricted signal clearly exceeds parallel and free."""
    ok, validated, msgs = True, 0, []
    for row in rows:
        if row["n_par"] == 0 or row["n_per"] == 0:
            msgs.append(f"  b={row['b']:.0f}: skipped (par n={row['n_par']}, perp n={row['n_per']})")
            continue
        validated += 1
        aniso = row["per"] - row["par"]
        above_free = row["per"] - row["free"]
        good = (aniso > 0.05) and (above_free > 0.03) and (row["par"] < row["per"])
        ok = ok and good
        msgs.append(f"  b={row['b']:.0f}: par<z>={row['par']:.4f} perp={row['per']:.4f} "
                    f"free={row['free']:.4f}  perp-par={aniso:+.4f} perp-free={above_free:+.4f} "
                    f"{'OK' if good else 'FAIL'}")
    if validated < 2:
        ok = False
        msgs.append(f"  FAIL: only {validated} shell(s) had both direction groups (need >=2)")
    return ok, msgs


def compare_analytic(r, nrmse_max):
    """(worst_shell_nrmse, rows) of sim S/S0 vs the van Gelderen curve.

    Prefers the frozen dmipy curve; if it has not been generated, falls back to
    the live pure-numpy GPA (``oracles.cylinder_vangelderen`` -- the same model),
    so the test RUNS (numpy-only) instead of skipping. Compares the full
    per-measurement pattern within each shell (captures the angular dependence),
    so NRMSE is normalised by the shell's range.
    """
    try:
        curve = oracles.load_analytic(ANALYTIC_NAME)
    except FileNotFoundError:
        curve = oracles.cylinder_vangelderen(r.sc, RADIUS, r.D)
    if curve.size != r.n:
        raise AssertionError(f"analytic curve ({curve.size}) != scheme rows ({r.n})")
    sc = r.sc
    nb = ~sc.b0_mask
    b_mm = np.round(sc.bvals / 1e6, 0)
    worst, rows = 0.0, []
    for b in sorted(set(b_mm[nb])):
        if b <= 0:
            continue
        m = nb & (b_mm == b)
        e = metrics.nrmse(r.ss0[m], curve[m], norm="range")
        worst = max(worst, e)
        rows.append((float(b), float(e)))
    return worst, rows


# --- pytest entry points (validation-analytic tier; needs pytest) ------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    return binary


def test_cylinder_anisotropy_cpu():
    """CPU: z-cylinder shows restricted-perpendicular vs free-parallel anisotropy."""
    ok, msgs = check_anisotropy(anisotropy_rows(run_cylinder(_binary_or_skip(), n=20000, gpu=0)))
    assert ok, "cylinder anisotropy failed:\n" + "\n".join(msgs)


def test_cylinder_anisotropy_gpu():
    """GPU: same anisotropy (self-skips without CUDA)."""
    import pytest
    r = run_cylinder(_binary_or_skip(), n=20000, gpu=1)
    if r.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    ok, msgs = check_anisotropy(anisotropy_rows(r))
    assert ok, "cylinder anisotropy failed (GPU):\n" + "\n".join(msgs)


def test_cylinder_vs_vangelderen():
    """CPU: signal matches the van Gelderen curve (frozen dmipy, else live GPA)."""
    r = run_cylinder(_binary_or_skip(), n=50000, t=1000, gpu=0)
    worst, _ = compare_analytic(r, nrmse_max=0.10)
    assert worst <= 0.10, f"worst per-shell NRMSE {worst:.4f} > 0.10 vs van Gelderen curve"


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--t", type=int, default=None)
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    ap.add_argument("--nrmse", type=float, default=0.10,
                    help="max allowed worst per-shell NRMSE vs the frozen analytic curve")
    a = ap.parse_args()

    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1
    r = run_cylinder(binary, n=a.n, t=a.t, gpu=a.gpu)
    if r.skipped:
        print("SKIP: binary built without CUDA / no GPU available.")
        return 0

    rows = anisotropy_rows(r)
    aniso_ok, msgs = check_anisotropy(rows)
    print(f"[cylinder anisotropy | N={a.n} {'GPU' if a.gpu else 'CPU'}]")
    for m in msgs:
        print(m)

    worst, arows = compare_analytic(r, a.nrmse)
    analytic_ok = worst <= a.nrmse
    print(f"[vs van Gelderen curve (frozen dmipy, else live GPA)] worst per-shell "
          f"NRMSE={worst:.4f} (max {a.nrmse}): {'OK' if analytic_ok else 'FAIL'}")
    for b, e in arows:
        print(f"  b={b:.0f}: NRMSE={e:.4f}")

    passed = aniso_ok and analytic_ok
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
