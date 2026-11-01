#!/usr/bin/env python3
"""
Cross-representation validation: the SAME cylinder, two geometry back-ends.

This is the dmipy-free cornerstone of the suite. A single radius-5 um z-cylinder
is simulated two ways and the diffusion signals are required to agree:

  (i)  NATIVE list  -- cylinder_restricted.conf: the simulator's analytic
       infinite Cylinder (Cylinder::checkCollision), walkers seeded intra.
  (ii) CAPPED PLY   -- cylinder_ply.conf: the SAME cylinder as a triangulated
       finite mesh (cylinder_capped.ply, r=5 um, caps at z=+/-200 um). Walkers
       are seeded in a central <spawning_area> (+/-50 um in z) so over the 50 ms
       diffusion time (1D RMS step ~7.7 um) they never reach the caps -> the mesh
       behaves like an infinite cylinder (pure radial restriction).

Because it is the SAME geometry with the SAME physics (N, T, D, scheme, seed),
the two S/S0 curves must coincide to within MC noise plus small representation
differences (a 128-facet polygon vs a smooth circle). The check is therefore
tight: per-shell mean |S/S0_list - S/S0_ply| within a narrow band. This
exercises the analytic-cylinder and the mesh collision paths against each other
with NO external oracle -- a wrong ply_scale, a wrong radius, or walkers reaching
the caps would all move S/S0 well beyond the band.

Two secondary checks:
  * HEX packing sanity (cylinder_hex_impermeable.conf): a DIFFERENT substrate
    (hex-packed cylinders, its own ICVF/topology), so we do NOT assert equality
    with the single cylinder -- only that it is physically restricted, i.e. the
    perpendicular S/S0 sits clearly above the free curve at multiple shells.
  * vs the frozen dmipy van Gelderen curve (skips cleanly if not present).

Tolerances are statistical (fixed-step MC noise ~1/sqrt(N) + O(1/T) bias); this
validates physics equivalence, not bit-exactness. Runs as pytest (CPU + a GPU
variant that self-skips without CUDA) OR standalone (numpy only, no pytest):

  python3 tests/validation/test_cross_representation.py --gpu 0 --n 8000
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

LIST_CONF = "cylinder_restricted.conf"   # native analytic infinite cylinder
PLY_CONF = "cylinder_ply.conf"           # capped triangulated mesh (same geometry)
HEX_CONF = "cylinder_hex_impermeable.conf"  # DIFFERENT substrate (restriction sanity)
ANALYTIC_NAME = "cylinder_r5um_d0p6"     # frozen van Gelderen curve (optional)
RADIUS = 5e-6                            # cylinder radius (matches both fixtures + ANALYTIC_NAME)
DIFFUSIVITY = 0.6e-9                      # m^2/s (cylinder_restricted.conf / cylinder_ply.conf)


def _vangelderen_curve(sc):
    """Frozen dmipy van Gelderen curve if present, else the live pure-numpy GPA.

    Mirrors test_cylinder_vangelderen.py so the PLY-vs-analytic cross-check RUNS
    dmipy-free (was: skip when the .npy is absent).
    """
    try:
        return oracles.load_analytic(ANALYTIC_NAME)
    except FileNotFoundError:
        return oracles.cylinder_vangelderen(sc, RADIUS, DIFFUSIVITY)

DIFFUSIVITY = 0.6e-9   # m^2/s -- identical in all three configs
PAR_MAX_DEG = 30.0     # angle-to-z <= 30 deg -> gradient ~along axis (free axial diffusion)
PERP_MIN_DEG = 60.0    # angle-to-z >= 60 deg -> gradient ~perpendicular (restricted)

# --- tolerances (measured deterministically at N=8000, seed 12345, num_process 1) --
# Observed worst per-shell mean |list-ply| = 0.0017, worst per-measurement = 0.0114.
# The bands below are generous multiples of that: they absorb MC noise (which grows
# ~1/sqrt(N) if N is lowered from the CLI) and the polygon-vs-circle discretization,
# while still catching any REAL representation error -- a wrong ply_scale/radius or
# walkers reaching the caps would shift S/S0 by >>0.02.
SHELL_TOL = 0.02       # per-shell mean |S/S0_list - S/S0_ply|
MEAS_TOL = 0.05        # per-measurement outlier guard
PERP_ABOVE_FREE = 0.03  # min perp(S/S0) - free for the "restricted" gate

Result = namedtuple("Result", "ss0 sc skipped")


def run_conf(binary, conf_name, n, gpu=0, t=None, timeout=600):
    """Run an accuracy config; return a :class:`Result` (normalised S/S0)."""
    conf_text = open(paths.accuracy(conf_name)).read()
    N = int(n)
    sc = scheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))
    overrides = {"N": N}
    if t:
        overrides["T"] = int(t)
    work = tempfile.mkdtemp(prefix="mcdc_xrep_")
    try:
        prefix = run.run_case(binary, conf_text, work, "x", paths.REPO_ROOT,
                              overrides=overrides, gpu=gpu, timeout=timeout)
        if prefix is None:
            return Result(None, sc, True)
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=N)
        if dwi.size != sc.n:
            raise AssertionError(f"{conf_name}: DWI entries ({dwi.size}) != scheme rows ({sc.n})")
        return Result(metrics.normalize_by_b0(dwi, sc.b0_mask), sc, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def agreement_rows(r_list, r_ply):
    """(rows, worst_meas): per-shell (b>0) mean S/S0 list vs ply + their abs diff.

    ``worst_meas`` is the largest per-measurement |list-ply| over the non-b0 rows
    (a loose outlier guard on top of the per-shell physical check).
    """
    sc = r_list.sc
    nb = ~sc.b0_mask
    ms_l = metrics.per_shell_mean(r_list.ss0, sc.bvals)
    ms_p = metrics.per_shell_mean(r_ply.ss0, sc.bvals)
    rows = [dict(b=b, lst=ms_l[b], ply=ms_p[b], d=abs(ms_l[b] - ms_p[b]))
            for b in sorted(k for k in ms_l if k > 0)]
    worst_meas = float(np.max(np.abs(r_list.ss0[nb] - r_ply.ss0[nb])))
    return rows, worst_meas


def check_agreement(rows, worst_meas, shell_tol=SHELL_TOL, meas_tol=MEAS_TOL):
    """(ok, messages): every shell mean agrees within ``shell_tol`` AND the worst
    per-measurement deviation stays within ``meas_tol``."""
    ok, msgs = True, []
    for row in rows:
        good = row["d"] <= shell_tol
        ok = ok and good
        msgs.append(f"  b={row['b']:.0f}: list={row['lst']:.4f} ply={row['ply']:.4f} "
                    f"|d|={row['d']:.4f} (tol {shell_tol}) {'OK' if good else 'FAIL'}")
    meas_ok = worst_meas <= meas_tol
    ok = ok and meas_ok
    msgs.append(f"  worst per-measurement |list-ply|={worst_meas:.4f} "
                f"(tol {meas_tol}) {'OK' if meas_ok else 'FAIL'}")
    return ok, msgs


def restriction_rows(r, D=DIFFUSIVITY):
    """Per-shell (b>0) perpendicular / parallel mean S/S0 + free exp(-bD)."""
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
                         free=math.exp(-(b * 1e6) * D)))
    return rows


def check_restricted(rows, above_free=PERP_ABOVE_FREE):
    """(ok, messages): perpendicular restricted signal clearly ABOVE free at >=2
    shells. A pure physics gate -- it makes NO claim that this substrate equals
    the single cylinder (the hex packing has its own radius/ICVF/topology)."""
    ok, validated, msgs = True, 0, []
    for row in rows:
        if row["n_per"] == 0:
            msgs.append(f"  b={row['b']:.0f}: skipped (perp n=0)")
            continue
        validated += 1
        gap = row["per"] - row["free"]
        good = gap > above_free
        ok = ok and good
        msgs.append(f"  b={row['b']:.0f}: perp={row['per']:.4f} par={row['par']:.4f} "
                    f"free={row['free']:.4f}  perp-free={gap:+.4f} {'OK' if good else 'FAIL'}")
    if validated < 2:
        ok = False
        msgs.append(f"  FAIL: only {validated} shell(s) had a perpendicular group (need >=2)")
    return ok, msgs


def compare_analytic(r, curve):
    """(worst_shell_nrmse, rows) of sim S/S0 vs a frozen van Gelderen ``curve``."""
    if curve.size != r.sc.n:
        raise AssertionError(f"analytic curve ({curve.size}) != scheme rows ({r.sc.n})")
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


def test_cross_representation_cpu():
    """CPU: native-list cylinder and capped-PLY cylinder give the same S/S0."""
    binary = _binary_or_skip()
    rl = run_conf(binary, LIST_CONF, n=8000, gpu=0)
    rp = run_conf(binary, PLY_CONF, n=8000, gpu=0)
    ok, msgs = check_agreement(*agreement_rows(rl, rp))
    assert ok, "list vs ply cross-representation mismatch:\n" + "\n".join(msgs)


def test_cross_representation_gpu():
    """GPU: same list-vs-ply equivalence (self-skips without CUDA)."""
    import pytest
    binary = _binary_or_skip()
    rl = run_conf(binary, LIST_CONF, n=8000, gpu=1)
    if rl.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    rp = run_conf(binary, PLY_CONF, n=8000, gpu=1)
    if rp.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    ok, msgs = check_agreement(*agreement_rows(rl, rp))
    assert ok, "list vs ply cross-representation mismatch (GPU):\n" + "\n".join(msgs)


def test_hex_restricted_cpu():
    """CPU: hex cylinder packing is physically restricted (perp S/S0 > free)."""
    r = run_conf(_binary_or_skip(), HEX_CONF, n=8000, gpu=0)
    ok, msgs = check_restricted(restriction_rows(r))
    assert ok, "hex packing not restricted:\n" + "\n".join(msgs)


def test_hex_restricted_gpu():
    """GPU: same hex restriction sanity (self-skips without CUDA)."""
    import pytest
    r = run_conf(_binary_or_skip(), HEX_CONF, n=8000, gpu=1)
    if r.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    ok, msgs = check_restricted(restriction_rows(r))
    assert ok, "hex packing not restricted (GPU):\n" + "\n".join(msgs)


def test_cross_vs_vangelderen():
    """CPU: BOTH the native-list and the capped-PLY cylinder match van Gelderen.

    The analytic curve is the frozen dmipy one if present, else the live pure-numpy
    GPA — so this PLY-vs-analytic check RUNS dmipy-free (the only test that ties a
    PLY mesh directly to an analytic formula).
    """
    binary = _binary_or_skip()
    rl = run_conf(binary, LIST_CONF, n=50000, t=1000, gpu=0)
    rp = run_conf(binary, PLY_CONF, n=50000, t=1000, gpu=0)
    curve = _vangelderen_curve(rl.sc)
    worst_l, _ = compare_analytic(rl, curve)
    worst_p, _ = compare_analytic(rp, curve)
    assert worst_l <= 0.10, f"list worst per-shell NRMSE {worst_l:.4f} > 0.10 vs van Gelderen"
    assert worst_p <= 0.10, f"ply worst per-shell NRMSE {worst_p:.4f} > 0.10 vs van Gelderen"


# --- standalone CLI (numpy only, no pytest) ----------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=8000, help="walkers per run (default 8000)")
    ap.add_argument("--t", type=int, default=None)
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    ap.add_argument("--shell-tol", type=float, default=SHELL_TOL)
    ap.add_argument("--meas-tol", type=float, default=MEAS_TOL)
    a = ap.parse_args()

    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1

    rl = run_conf(binary, LIST_CONF, n=a.n, t=a.t, gpu=a.gpu)
    if rl.skipped:
        print("SKIP: binary built without CUDA / no GPU available.")
        return 0
    rp = run_conf(binary, PLY_CONF, n=a.n, t=a.t, gpu=a.gpu)
    if rp.skipped:
        print("SKIP: binary built without CUDA / no GPU available.")
        return 0

    rows, worst_meas = agreement_rows(rl, rp)
    agree_ok, msgs = check_agreement(rows, worst_meas, a.shell_tol, a.meas_tol)
    print(f"[cross-representation: list vs ply | N={a.n} {'GPU' if a.gpu else 'CPU'}]")
    for m in msgs:
        print(m)

    rhex = run_conf(binary, HEX_CONF, n=a.n, t=a.t, gpu=a.gpu)
    hex_ok, hmsgs = check_restricted(restriction_rows(rhex))
    print(f"[hex packing restriction sanity | N={a.n}]")
    for m in hmsgs:
        print(m)

    curve = _vangelderen_curve(rl.sc)
    worst_l, _ = compare_analytic(rl, curve)
    worst_p, _ = compare_analytic(rp, curve)
    van_ok = worst_l <= 0.10 and worst_p <= 0.10
    print(f"[vs van Gelderen (frozen dmipy, else live GPA)] list NRMSE={worst_l:.4f} "
          f"ply NRMSE={worst_p:.4f} (max 0.10): {'OK' if van_ok else 'FAIL'}")

    passed = agree_ok and hex_ok and van_ok
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
