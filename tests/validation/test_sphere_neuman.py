#!/usr/bin/env python3
"""
Analytical validation: restricted diffusion in a single impermeable sphere.

The SPH-chain analogue of test_cylinder_vangelderen.py. It runs
tests/accuracy/sphere_restricted.conf (one native-list sphere, radius 5 um at the
origin, walkers seeded intra, IMPERMEABLE, T2 ~ inf) and validates the restricted
signal two ways:

(a) RESTRICTION + ISOTROPY SANITY (always on, no dmipy needed).
    A sphere restricts diffusion in EVERY direction, so at every shell the mean
    S/S0 must sit clearly ABOVE the free curve exp(-b*D) (restriction raises the
    signal), and -- because the sphere is isotropic -- the per-direction scatter
    within a shell must be tiny (all directions share one value). These two cheap
    checks (raised + isotropic) together are a robust physics gate that needs no
    external oracle.

(b) vs the frozen Neuman / Murday-Cotts curve (skips if the curve is absent).
    Compares the per-shell signal to the dmipy S4SphereGaussianPhaseApproximation
    curve frozen at tests/accuracy/analytic/ (see tests/tools/gen_analytic.py
    --model sphere). Because the sphere is isotropic the in-shell range is ~0, so
    per-shell NRMSE is normalised by the shell MEAN (not its range, which would be
    degenerate). Valid only at small b (the GPA is a closed impermeable sphere).

Tolerances are statistical: fixed-step MC has noise ~1/sqrt(N) plus an O(1/T)
finite-step bias that grows with b -- this validates the physics, not
bit-exactness. Runs as pytest (CPU + GPU-self-skip) or standalone:

  python3 tests/validation/test_sphere_neuman.py --gpu 0 --n 20000
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

CONF = paths.accuracy("sphere_restricted.conf")
ANALYTIC_NAME = "sphere_r5um_d0p6"
RADIUS = 5e-6         # single_sphere.list radius (matches ANALYTIC_NAME r5um)

# Physics gate thresholds (see explore/verify comments):
#   restricted mean must exceed free by this absolute margin AND ratio, on >=2 shells;
#   observed on the current build: margins 0.18..0.41, ratios 1.24..2.37 -> huge headroom.
ABOVE_FREE_MIN = 0.03     # <S/S0>_restricted - exp(-bD)  must exceed this
RATIO_FREE_MIN = 1.05     # <S/S0>_restricted / exp(-bD)  must exceed this
ISO_CV_MAX = 0.05         # per-shell std/mean (isotropy): observed ~0.001..0.002

Result = namedtuple("Result", "ss0 sc D n skipped")


def run_sphere(binary, n=20000, t=None, gpu=0, timeout=600):
    """Run the single-sphere case; return a :class:`Result` (normalised S/S0)."""
    conf_text = open(CONF).read()
    N = int(n) if n else int(io.read_conf_value(conf_text, "N"))
    D = float(io.read_conf_value(conf_text, "diffusivity"))
    sc = scheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))
    overrides = {"N": N}
    if t:
        overrides["T"] = int(t)
    work = tempfile.mkdtemp(prefix="mcdc_sph_")
    try:
        prefix = run.run_case(binary, conf_text, work, "sph", paths.REPO_ROOT,
                              overrides=overrides, gpu=gpu, timeout=timeout)
        if prefix is None:
            return Result(None, None, None, 0, True)
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=N)
        if dwi.size != sc.n:
            raise AssertionError(f"DWI entries ({dwi.size}) != scheme rows ({sc.n})")
        return Result(metrics.normalize_by_b0(dwi, sc.b0_mask), sc, D, sc.n, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def restriction_rows(r):
    """Per-shell (b>0): mean S/S0, its std (isotropy), and the free exp(-bD)."""
    sc = r.sc
    nb = ~sc.b0_mask
    b_mm = np.round(sc.bvals / 1e6, 0)
    rows = []
    for b in sorted(set(b_mm[nb])):
        if b <= 0:
            continue
        m = nb & (b_mm == b)
        vals = r.ss0[m]
        mean = float(vals.mean())
        rows.append(dict(b=float(b), n=int(vals.size), mean=mean, std=float(vals.std()),
                         free=math.exp(-(b * 1e6) * r.D)))
    return rows


def check_restriction(rows):
    """(ok, messages). Require >=2 shells where the restricted mean clearly
    exceeds the free curve AND the in-shell scatter is tiny (isotropy)."""
    ok, validated, msgs = True, 0, []
    for row in rows:
        if row["n"] == 0:
            continue
        validated += 1
        above = row["mean"] - row["free"]
        ratio = row["mean"] / row["free"]
        cv = row["std"] / row["mean"] if row["mean"] > 0 else float("inf")
        good = (above > ABOVE_FREE_MIN) and (ratio > RATIO_FREE_MIN) and (cv < ISO_CV_MAX)
        ok = ok and good
        msgs.append(f"  b={row['b']:.0f}: <S/S0>={row['mean']:.4f} free={row['free']:.4f} "
                    f"above={above:+.4f} ratio={ratio:.3f} iso_cv={cv:.4f} "
                    f"{'OK' if good else 'FAIL'}")
    if validated < 2:
        ok = False
        msgs.append(f"  FAIL: only {validated} shell(s) validated (need >=2)")
    return ok, msgs


def compare_analytic(r, nrmse_max):
    """(worst_shell_nrmse, rows) of sim S/S0 vs the Neuman sphere curve.

    Prefers the frozen dmipy curve; if it has not been generated, falls back to
    the live pure-numpy GPA (``oracles.sphere_neuman`` -- the same model), so the
    test RUNS (numpy-only) instead of skipping. The sphere is isotropic (in-shell
    range ~0), so NRMSE is normalised by the shell MEAN -- range-normalisation
    would divide by ~0 and blow up.
    """
    try:
        curve = oracles.load_analytic(ANALYTIC_NAME)
    except FileNotFoundError:
        curve = oracles.sphere_neuman(r.sc, RADIUS, r.D)
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
        e = metrics.nrmse(r.ss0[m], curve[m], norm="mean")
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


def test_sphere_restriction_cpu():
    """CPU: sphere signal is isotropic and sits above free exp(-bD) (restriction)."""
    ok, msgs = check_restriction(restriction_rows(run_sphere(_binary_or_skip(), n=20000, gpu=0)))
    assert ok, "sphere restriction/isotropy failed:\n" + "\n".join(msgs)


def test_sphere_restriction_gpu():
    """GPU: same restriction/isotropy (self-skips without CUDA)."""
    import pytest
    r = run_sphere(_binary_or_skip(), n=20000, gpu=1)
    if r.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    ok, msgs = check_restriction(restriction_rows(r))
    assert ok, "sphere restriction/isotropy failed (GPU):\n" + "\n".join(msgs)


def test_sphere_vs_neuman():
    """CPU: signal matches the Neuman sphere curve (frozen dmipy, else live GPA)."""
    r = run_sphere(_binary_or_skip(), n=50000, t=1000, gpu=0)
    worst, _ = compare_analytic(r, nrmse_max=0.10)
    assert worst <= 0.10, f"worst per-shell NRMSE {worst:.4f} > 0.10 vs Neuman sphere curve"


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
    r = run_sphere(binary, n=a.n, t=a.t, gpu=a.gpu)
    if r.skipped:
        print("SKIP: binary built without CUDA / no GPU available.")
        return 0

    rows = restriction_rows(r)
    sane_ok, msgs = check_restriction(rows)
    print(f"[sphere restriction + isotropy | N={a.n} {'GPU' if a.gpu else 'CPU'}]")
    for m in msgs:
        print(m)

    worst, arows = compare_analytic(r, a.nrmse)
    analytic_ok = worst <= a.nrmse
    print(f"[vs Neuman sphere curve (frozen dmipy, else live GPA)] worst per-shell "
          f"NRMSE={worst:.4f} (max {a.nrmse}): {'OK' if analytic_ok else 'FAIL'}")
    for b, e in arows:
        print(f"  b={b:.0f}: NRMSE={e:.4f}")

    passed = sane_ok and analytic_ok
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
