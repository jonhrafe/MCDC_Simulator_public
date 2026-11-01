#!/usr/bin/env python3
"""
Validation: two-compartment PERMEABLE exchange + the Kärger oracle.

Three parts, of increasing dependence:

(1) ORACLE SELF-TEST (always on, numpy only). oracles.karger is the closed-form
    Kärger two-site exchange signal. We check it against its defining physics: the
    b=0 value, the NO-exchange limit (static bi-exponential), the FAST-exchange
    limit (mono-exponential at the mean D), and -- the strong check -- a direct
    matrix exponential of the 2x2 Kärger generator computed independently with
    numpy. This validates the oracle regardless of the simulator.

(2) FIXTURE SANITY (always on, needs the binary). Runs
    tests/accuracy/karger_exchange.conf (one permeable sphere, distinct intra/extra
    D, slow exchange, walkers seeded uniformly) and checks it is a well-formed
    two-compartment signal: S/S0 in (0,1], monotonically decaying in b, and a b0
    intra fraction matching the KNOWN geometric fraction V_sphere/V_voxel.

(3) QUANTITATIVE KÄRGER FIT (SKIPPED, with a TODO). A raw bi-exponential Kärger
    fit of this MC fixture is deliberately NOT asserted: intra diffusion in a 5 um
    sphere is RESTRICTED (non-Gaussian), which violates Kärger's Gaussian-
    compartment assumption. See the skip message for the geometry/regime a real
    quantitative Kärger validation would require. (Per the suite policy: prefer a
    correct skip over a passing-but-physically-wrong assertion.)

  python3 tests/validation/test_karger_exchange.py --n 8000 --gpu 0
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

CONF = paths.accuracy("karger_exchange.conf")

SKIP_TODO = (
    "quantitative Kärger fit not asserted: intra diffusion in a 5 um sphere is "
    "RESTRICTED (non-Gaussian), so the apparent intra D is time/b-dependent, not "
    "d_intra -- Kärger assumes Gaussian compartments. A real Kärger validation "
    "needs EITHER (i) Gaussian compartments (large/weakly-restricting compartments, "
    "or periodic slabs probed ALONG the membrane-parallel free direction) OR "
    "(ii) apparent per-compartment D fitted from the impermeable (kappa=0) runs "
    "and fed to oracles.karger, in the slow/intermediate-exchange regime "
    "(per-encounter crossing prob << 1, exchange observable within the diffusion "
    "time). Until then: oracle validated by self-test, fixture by geometry + sanity."
)


# --------------------------- oracle self-test --------------------------------
def _matrix_karger(bvals, f_a, D_a, D_b, k):
    """Independent reference: E(b) = [1,1] . expm(-b*L) . [f_a, f_b] with the
    Kärger generator L = [[Da+ka/b, -kb/b], [-ka/b, Db+kb/b]] (ka=k*fb, kb=k*fa),
    via a numpy eigendecomposition -- a path fully separate from oracles.karger's
    eigenvalue/weight algebra, so agreement cross-validates that algebra."""
    f_b = 1.0 - f_a
    k_a, k_b = k * f_b, k * f_a
    out = []
    for b in np.asarray(bvals, dtype=np.float64):
        if b <= 0.0:
            out.append(1.0)
            continue
        L = np.array([[D_a + k_a / b, -k_b / b],
                      [-k_a / b, D_b + k_b / b]], dtype=np.float64)
        w, V = np.linalg.eig(-b * L)
        expm = (V * np.exp(w)) @ np.linalg.inv(V)
        out.append(float(np.real(np.array([1.0, 1.0]) @ expm @ np.array([f_a, f_b]))))
    return np.array(out, dtype=np.float64)


def karger_selftest():
    """(ok, messages) -- validate oracles.karger against limits + matrix expm."""
    b = np.array([0.0, 500.0, 1000.0, 2000.0, 4000.0]) * 1e6      # s/m^2
    f_a, D_a, D_b = 0.3, 0.6e-9, 2.0e-9
    f_b, Dbar = 1.0 - f_a, f_a * 0.6e-9 + 0.7 * 2.0e-9
    msgs, ok = [], True

    e0 = oracles.karger(b, f_a, D_a, D_b, k=20.0)
    d0 = abs(e0[0] - 1.0)
    ok &= d0 < 1e-12
    msgs.append(f"  S(b=0)=1: dev={d0:.2e} {'OK' if d0 < 1e-12 else 'FAIL'}")

    no_ex = oracles.karger(b, f_a, D_a, D_b, k=0.0)
    ref_no = f_a * np.exp(-b * D_a) + f_b * np.exp(-b * D_b)
    dno = float(np.max(np.abs(no_ex - ref_no)))
    ok &= dno < 1e-12
    msgs.append(f"  no-exchange -> bi-exp: max dev={dno:.2e} {'OK' if dno < 1e-12 else 'FAIL'}")

    fast = oracles.karger(b, f_a, D_a, D_b, k=1e7)
    ref_fast = np.exp(-b * Dbar)
    dfast = float(np.max(np.abs(fast - ref_fast)))
    ok &= dfast < 5e-4
    msgs.append(f"  fast-exchange -> mono-exp(Dbar): max dev={dfast:.2e} "
                f"{'OK' if dfast < 5e-4 else 'FAIL'}")

    for k in (5.0, 20.0, 100.0):
        e = oracles.karger(b, f_a, D_a, D_b, k=k)
        em = _matrix_karger(b, f_a, D_a, D_b, k=k)
        dm = float(np.max(np.abs(e - em)))
        ok &= dm < 1e-9
        msgs.append(f"  vs matrix expm (k={k:g}/s): max dev={dm:.2e} {'OK' if dm < 1e-9 else 'FAIL'}")
    return ok, msgs


# --------------------------- fixture sanity ----------------------------------
Result = namedtuple("Result", "ss0 sc intra_frac geom_frac n skipped")


def _geom_intra_fraction():
    """Known f_intra = V_sphere / V_voxel from the conf voxel + the list radius."""
    text = open(CONF).read()
    lines = text.splitlines()
    vi = next(i for i, ln in enumerate(lines) if ln.strip() == "<voxels>")
    lo = [float(x) for x in lines[vi + 1].split()]
    hi = [float(x) for x in lines[vi + 2].split()]
    v_vox = abs((hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]))
    lst = os.path.join(paths.REPO_ROOT, io.read_conf_value(text, "spheres_list"))
    rows = [ln for ln in open(lst).read().splitlines()
            if ln.strip() and not ln.strip().startswith("#")]
    scale = float(rows[0].split()[0])
    radius = float(rows[1].split()[3]) * scale        # `x y z r p t2`, r column, * scale
    v_sph = 4.0 / 3.0 * math.pi * radius ** 3
    return v_sph / v_vox


def run_fixture(binary, n=None, t=None, gpu=0, timeout=600):
    conf_text = open(CONF).read()
    N = int(n) if n else int(io.read_conf_value(conf_text, "N"))
    sc = scheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))
    overrides = {"N": N}
    if t:
        overrides["T"] = int(t)
    work = tempfile.mkdtemp(prefix="mcdc_karger_")
    try:
        prefix = run.run_case(binary, conf_text, work, "karger", paths.REPO_ROOT,
                              overrides=overrides, gpu=gpu, timeout=timeout)
        if prefix is None:
            return Result(None, None, None, None, 0, True)
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=N)
        intra = io.load_dwi(prefix + "_DWI_intra.txt", n_walkers=N)
        if dwi.size != sc.n:
            raise AssertionError(f"DWI entries ({dwi.size}) != scheme rows ({sc.n})")
        ss0 = metrics.normalize_by_b0(dwi, sc.b0_mask)
        intra_frac = float(intra[sc.b0_mask].mean())     # T2 inf -> b0 == final-location count/N
        return Result(ss0, sc, intra_frac, _geom_intra_fraction(), sc.n, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def check_fixture(r):
    """(ok, messages) -- well-formed decaying two-compartment signal, known frac."""
    sc = r.sc
    nb = ~sc.b0_mask
    b_mm = np.round(sc.bvals / 1e6, 0)
    means = [float(r.ss0[nb & (b_mm == b)].mean())
             for b in sorted(set(b_mm[nb])) if b > 0]
    msgs, ok = [], True

    finite = bool(np.all(np.isfinite(r.ss0)))
    bounded = bool(np.all(r.ss0 > 0) and np.all(r.ss0 <= 1.0 + 1e-3))
    monotone = all(x > y for x, y in zip(means, means[1:]))
    frac_ok = abs(r.intra_frac - r.geom_frac) < 0.05
    ok = finite and bounded and monotone and frac_ok

    msgs.append(f"  per-shell mean S/S0: {[round(m, 4) for m in means]}")
    msgs.append(f"  finite={finite} bounded(0,1]={bounded} monotone_decay={monotone}")
    msgs.append(f"  intra fraction: sim={r.intra_frac:.4f} geom={r.geom_frac:.4f} "
                f"(|d|<0.05: {frac_ok})")
    return ok, msgs


# --- pytest entry points -----------------------------------------------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    return binary


def test_karger_oracle_selftest():
    """Always-on, numpy-only: oracles.karger matches its limits + matrix expm."""
    ok, msgs = karger_selftest()
    assert ok, "karger oracle self-test failed:\n" + "\n".join(msgs)


def test_karger_fixture_sanity_cpu():
    """CPU: the permeable two-compartment fixture is a valid decaying signal
    with the known geometric intra fraction."""
    ok, msgs = check_fixture(run_fixture(_binary_or_skip(), n=8000, gpu=0))
    assert ok, "karger fixture sanity failed:\n" + "\n".join(msgs)


def test_karger_fixture_sanity_gpu():
    """GPU: same fixture sanity (self-skips without CUDA)."""
    import pytest
    r = run_fixture(_binary_or_skip(), n=8000, gpu=1)
    if r.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    ok, msgs = check_fixture(r)
    assert ok, "karger fixture sanity failed (GPU):\n" + "\n".join(msgs)


def test_karger_quantitative_fit():
    """SKIP: quantitative bi-exponential Kärger fit -- not valid for restricted
    intra diffusion. See SKIP_TODO for the geometry/regime required."""
    import pytest
    pytest.skip(SKIP_TODO)


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=8000)
    ap.add_argument("--t", type=int, default=None)
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    a = ap.parse_args()

    print("[karger oracle self-test | numpy only]")
    oracle_ok, msgs = karger_selftest()
    for m in msgs:
        print(m)
    print(f"  -> {'PASS' if oracle_ok else 'FAIL'}")

    binary = paths.binary_path(a.bin)
    fixture_ok = True
    if not os.path.exists(binary):
        print(f"[fixture sanity] SKIP: simulator binary not found: {binary}")
    else:
        r = run_fixture(binary, n=a.n, t=a.t, gpu=a.gpu)
        if r.skipped:
            print("[fixture sanity] SKIP: binary built without CUDA / no GPU available.")
        else:
            fixture_ok, fmsgs = check_fixture(r)
            print(f"[fixture sanity | N={a.n} {'GPU' if a.gpu else 'CPU'}]")
            for m in fmsgs:
                print(m)
            print(f"  -> {'PASS' if fixture_ok else 'FAIL'}")

    print("[quantitative Kärger fit] SKIP: " + SKIP_TODO.split(":", 1)[0])
    ok = oracle_ok and fixture_ok
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
