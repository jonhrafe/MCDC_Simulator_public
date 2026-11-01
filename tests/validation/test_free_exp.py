#!/usr/bin/env python3
"""
Analytical validation: free diffusion  ->  S/S0 = exp(-b*D).

The simplest oracle in the suite and the template for every other analytical
test: run a fixed config, parse the scheme, compute the closed-form signal per
measurement, and assert agreement within a statistical band.

With NO obstacles the PGSE signal must follow the free-diffusion law
``S/S0 = exp(-b*D)`` for every measurement (b from the scheme, D from the conf,
both SI). The b0 signal equals N (all phases zero, T2=inf), so ``S/N`` is already
S/S0 -- we compare it directly to ``exp(-b*D)``.

Tolerance: fixed-length MC steps approximate Gaussian diffusion only to O(1/T),
a bias that grows with b (~1.7% at b=2000, T=1000), on top of MC noise ~1/sqrt(N).
So the band is a few % -- this validates the physics, not bit-exactness.

Runs as a pytest test (CPU + a GPU variant that self-skips without CUDA) OR as a
standalone script (used by CTest for the fast smoke check):

  python3 tests/validation/test_free_exp.py --n 500 --t 300 --gpu 0 --rtol 0.15
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

CONF = paths.accuracy("free_diffusion.conf")

# Free diffusion is isotropic, so the physical quantity is the PER-SHELL MEAN
# (all directions in a shell share exp(-b*D)); per-measurement scatter is pure MC
# noise. We validate the shell means and keep a loose per-measurement bound only
# as a gross-outlier guard.
Result = namedtuple("Result", "worst_meas worst_shell n skipped")


def run_free_exp(binary, n=None, t=None, gpu=0, timeout=600):
    """Run free diffusion; return a :class:`Result`.

    ``worst_shell`` = largest ``|<S/N> - exp(-bD)| / exp(-bD)`` over non-b0 shell
    means (the physical check). ``worst_meas`` = largest per-measurement relative
    deviation (a loose outlier guard). ``skipped`` is True iff the GPU path bailed.
    Raises on a real simulator error or a shape mismatch.
    """
    conf_text = open(CONF).read()
    N = int(n) if n else int(io.read_conf_value(conf_text, "N"))
    D = float(io.read_conf_value(conf_text, "diffusivity"))       # m^2/s
    sc = scheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))

    overrides = {"N": N}
    if t:
        overrides["T"] = int(t)
    work = tempfile.mkdtemp(prefix="mcdc_freeexp_")
    try:
        prefix = run.run_case(binary, conf_text, work, "free", paths.REPO_ROOT,
                              overrides=overrides, gpu=gpu, timeout=timeout)
        if prefix is None:
            return Result(None, None, 0, True)
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=N)       # S/N (= S/S0 here)
        if dwi.size != sc.n:
            raise AssertionError(f"DWI entries ({dwi.size}) != scheme rows ({sc.n})")
        expected = oracles.free_exp(sc.bvals, D)
        worst_meas = float(np.max(np.abs(dwi - expected) / expected))
        run_shell = metrics.per_shell_mean(dwi, sc.bvals)
        exp_shell = metrics.per_shell_mean(expected, sc.bvals)
        worst_shell = max(abs(run_shell[b] - exp_shell[b]) / exp_shell[b]
                          for b in exp_shell if b > 0.0)
        return Result(worst_meas, float(worst_shell), sc.n, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


# --- pytest entry points (validation-analytic tier; needs pytest) ------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    return binary


def _assert_free_exp(r):
    assert r.worst_shell is not None and not math.isnan(r.worst_shell), "NaN signal"
    assert r.worst_shell <= 0.03, f"per-shell mean dev {r.worst_shell:.4f} > 0.03 ({r.n} meas)"
    assert r.worst_meas <= 0.10, f"per-measurement outlier {r.worst_meas:.4f} > 0.10"


def test_free_exp_cpu():
    """CPU free diffusion tracks exp(-bD) within a statistical band."""
    _assert_free_exp(run_free_exp(_binary_or_skip(), n=50000, t=1000, gpu=0))


def test_free_exp_gpu():
    """GPU free diffusion tracks exp(-bD) (skips cleanly without CUDA)."""
    import pytest
    r = run_free_exp(_binary_or_skip(), gpu=1)
    if r.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    _assert_free_exp(r)


# --- diffusivity levels: recover ADC == D at low AND high D -------------------
# Meaningful levels (not a sweep): the step length is sqrt(6 D dt), so the recovered
# free-diffusion ADC must equal D whether D is small or large. Checked at the LOWEST
# b-shell (robust; a raw exp(-bD) match gets noisy once b*D is large).
DIFFUSIVITIES = (0.5e-9, 3.0e-9)      # low / high (m^2/s)


def free_adc(binary, D, n=50000, t=1000, gpu=0, timeout=600):
    """Run free diffusion at diffusivity ``D``; return recovered ADC at the lowest b (or None)."""
    conf_text = open(CONF).read()
    sc = scheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))
    work = tempfile.mkdtemp(prefix="mcdc_freeD_")
    try:
        prefix = run.run_case(binary, conf_text, work, "fd", paths.REPO_ROOT,
                              overrides={"N": n, "T": t, "diffusivity": repr(float(D))},
                              gpu=gpu, timeout=timeout)
        if prefix is None:
            return None
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=n)
        shell = metrics.per_shell_mean(dwi, sc.bvals)
        b = min(x for x in shell if x > 0.0)                 # lowest non-b0 shell (s/mm^2)
        return float(-math.log(max(shell[b], 1e-12)) / (b * 1e6))
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _assert_diffusivity(binary, gpu):
    import pytest
    for D in DIFFUSIVITIES:
        adc = free_adc(binary, D, gpu=gpu)
        if adc is None:
            pytest.skip("binary built without CUDA / no GPU available")
        assert abs(adc - D) / D <= 0.05, \
            f"free ADC {adc:.3e} != D {D:.3e} ({'GPU' if gpu else 'CPU'})"


def test_free_exp_diffusivity_cpu():
    """CPU: recovered free ADC equals D at both low and high diffusivity."""
    _assert_diffusivity(_binary_or_skip(), gpu=0)


def test_free_exp_diffusivity_gpu():
    """GPU: same (self-skips without CUDA)."""
    _assert_diffusivity(_binary_or_skip(), gpu=1)


# --- standalone CLI (smoke tier; numpy-only, no pytest) ----------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=None, help="override N (default: conf value)")
    ap.add_argument("--t", type=int, default=None, help="override T (default: conf value)")
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    ap.add_argument("--rtol", type=float, default=0.05,
                    help="max allowed per-shell-mean relative deviation")
    a = ap.parse_args()

    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1
    r = run_free_exp(binary, n=a.n, t=a.t, gpu=a.gpu)
    if r.skipped:
        print("SKIP: binary built without CUDA / no GPU available.")
        return 0
    if r.worst_shell is None or math.isnan(r.worst_shell):
        print("FAIL: NaN signal (no valid DWI).")
        return 1
    ok = r.worst_shell <= a.rtol
    print(f"{'PASS' if ok else 'FAIL'}: free diffusion vs exp(-bD): "
          f"per-shell mean dev={r.worst_shell:.4f} (rtol {a.rtol}), "
          f"per-meas worst={r.worst_meas:.4f} "
          f"[N={a.n or 'conf'}, {r.n} meas, {'GPU' if a.gpu else 'CPU'}]")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
