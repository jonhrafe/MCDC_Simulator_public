#!/usr/bin/env python3
"""
Analytical validation: T2 relaxation weighting -- CLOSED FORM, no dmipy.

The T2-weighted signal carries a factor exp(-TE/T2) at each measurement's echo
time TE, on top of the diffusion attenuation. The simulator accumulates this per
step (walker.t2_log[t] = t2_log[t-1]*exp(-dt/T2_eff)) and samples it at the echo,
using the compartment T2 the walker is in (intra: sphere T2; extra: t2_extra).
Crucially the T2 weighting is applied POST-HOC and never feeds back into the
walk, so with a fixed seed the trajectories -- and hence the intra/extra split --
are identical whether or not T2 decays. That makes a same-seed differential
exact and MC-noise-free.

Fixture: tests/accuracy/sphere_multiT2.conf (impermeable sphere, separate_signals
so we get the intra and extra signals separately; intra T2 from the list, extra
T2 = t2_extra). The scheme (PGSE_sample_scheme.scheme) has a SINGLE TE, so we
validate the b0 signal MAGNITUDE:

    run A: as-is                                     -> S_c(b0) = N_c exp(-TE/T2_c)
    run B: reference with T2 -> inf in BOTH comparts -> S_c(b0) = N_c
    =>  S_c^A(b0) / S_c^B(b0) = exp(-TE / T2_c)      for c in {intra, extra}

Both compartment ratios must equal exp(-TE/T2_c) for the KNOWN T2 (t2_extra from
the conf, in seconds; the sphere-list T2 column is also in seconds since v2.3.1,
like every other file quantity of an SI conf). We also check the
separate-signal identity S_tot(b0) = S_intra(b0) + S_extra(b0). (If a scheme with
MULTIPLE TEs is used, we instead fit ln S_extra(b0) vs TE -> slope = -1/T2.)

Tolerance: the ratio is deterministic (every same-compartment spin at b0 carries
the identical factor exp(-TE/T2), so there is no 1/sqrt(N) noise); the ~2% band
only covers the text-output precision, the int(TE/dt) echo-index flooring and the
odd boundary walker. Observed on the current build: extra ~1e-5, intra ~1e-3.

  python3 tests/validation/test_t2_decay.py --n 6000 --gpu 0 --rtol 0.02
"""
import argparse
import math
import os
import shutil
import sys
import tempfile
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcdc_testkit import io, paths, run, scheme  # noqa: E402

import numpy as np  # noqa: E402

CONF = paths.accuracy("sphere_multiT2.conf")
INF_T2 = "1000000000"     # >1e5 -> simulator treats the compartment T2 as ~infinite

Result = namedtuple("Result", "sc te t2_extra t2_intra "
                              "iA eA tA iB eB N multi skipped")


def _read_list_t2_seconds(spheres_list_rel):
    """First sphere's T2 from a spheres_list, in SECONDS.

    Since v2.3.1 the spheres_list T2 column is read in SECONDS in an SI conf (it used to
    be read raw, as internal ms), the same unit as the conf's t2_extra, so it is compared
    with TE (seconds) as is. Row layout: `x y z r p t2`.
    """
    src = os.path.join(paths.REPO_ROOT, spheres_list_rel)
    rows = [ln for ln in open(src).read().splitlines()
            if ln.strip() and not ln.strip().startswith("#")]
    return float(rows[1].split()[-1])             # rows[0] is the scale line


def _inf_t2_listfile(spheres_list_rel, workdir):
    """Write a copy of the spheres_list with every T2 set to ~infinity."""
    src = os.path.join(paths.REPO_ROOT, spheres_list_rel)
    out = []
    for i, ln in enumerate(open(src).read().splitlines()):
        s = ln.strip()
        if i == 0 or not s or s.startswith("#"):
            out.append(ln)
            continue
        toks = s.split()
        toks[-1] = INF_T2                          # T2 is the last column
        out.append(" ".join(toks))
    dst = os.path.join(workdir, "inf_t2.list")
    with open(dst, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return dst


def run_t2(binary, n=None, gpu=0, t2_extra_s=None, timeout=600):
    """Run the T2 fixture (A) and its T2->inf reference (B); return a Result.

    ``t2_extra_s`` overrides the extra-compartment T2 (seconds) -- pass distinct
    values (e.g. 0.040 and 0.150) to check exp(-TE/T2) at a low AND a high T2.
    """
    conf_text = open(CONF).read()
    N = int(n) if n else int(io.read_conf_value(conf_text, "N"))
    sc = scheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))
    t2_extra = float(t2_extra_s) if t2_extra_s is not None else float(io.read_conf_value(conf_text, "t2_extra"))
    spheres_list_rel = io.read_conf_value(conf_text, "spheres_list")
    t2_intra = _read_list_t2_seconds(spheres_list_rel)                       # seconds
    te_b0 = np.unique(np.round(sc.TE[sc.b0_mask], 6))
    multi = te_b0.size > 1
    a_over = {"N": N}
    if t2_extra_s is not None:
        a_over["t2_extra"] = repr(float(t2_extra_s))

    work = tempfile.mkdtemp(prefix="mcdc_t2_")
    try:
        pa = run.run_case(binary, conf_text, work, "t2a", paths.REPO_ROOT,
                          overrides=a_over, gpu=gpu, timeout=timeout)
        if pa is None:
            return Result(sc, None, t2_extra, t2_intra, None, None, None,
                          None, None, N, multi, True)
        iA = io.load_dwi(pa + "_DWI_intra.txt", n_walkers=N)
        eA = io.load_dwi(pa + "_DWI_extra.txt", n_walkers=N)
        tA = io.load_dwi(pa + "_DWI.txt", n_walkers=N)
        iB = eB = None
        if not multi:
            inf_list = _inf_t2_listfile(spheres_list_rel, work)
            pb = run.run_case(binary, conf_text, work, "t2b", paths.REPO_ROOT,
                              overrides={"N": N, "t2_extra": INF_T2, "spheres_list": inf_list},
                              gpu=gpu, timeout=timeout)
            iB = io.load_dwi(pb + "_DWI_intra.txt", n_walkers=N)
            eB = io.load_dwi(pb + "_DWI_extra.txt", n_walkers=N)
        te = float(te_b0[0]) if not multi else None
        return Result(sc, te, t2_extra, t2_intra, iA, eA, tA, iB, eB, N, multi, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def evaluate(r):
    """Return (ok, worst_rel_dev, messages) for the T2 checks."""
    b0 = r.sc.b0_mask
    msgs, worst = [], 0.0

    # separate-signal identity: S_tot(b0) == S_intra(b0) + S_extra(b0)
    tot, split = float(r.tA[b0].mean()), float((r.iA + r.eA)[b0].mean())
    cons = abs(tot - split) / max(tot, 1e-30)
    worst = max(worst, cons)
    msgs.append(f"  separate-signal identity: tot={tot:.5f} intra+extra={split:.5f} "
                f"rel={cons:.2e}")

    if r.multi:
        # multi-TE: fit ln S_extra(b0) vs TE -> slope = -1/T2_extra (removes N_extra).
        te = r.sc.TE[b0]
        s = r.eA[b0]
        te_u = np.unique(np.round(te, 6))
        y = np.array([np.log(s[np.round(te, 6) == t].mean()) for t in te_u])
        slope = np.polyfit(te_u, y, 1)[0]
        t2_fit = -1.0 / slope
        dev = abs(t2_fit - r.t2_extra) / r.t2_extra
        worst = max(worst, dev)
        msgs.append(f"  multi-TE extra fit: T2={t2_fit:.4f}s expected={r.t2_extra:.4f}s dev={dev:.4f}")
        return worst <= 0.05, worst, msgs

    # single-TE differential: S_c^A(b0)/S_c^B(b0) = exp(-TE/T2_c)
    for name, a, b, t2 in (("extra", r.eA, r.eB, r.t2_extra),
                           ("intra", r.iA, r.iB, r.t2_intra)):
        num, den = float(a[b0].mean()), float(b[b0].mean())
        ratio = num / den
        expected = math.exp(-r.te / t2)
        dev = abs(ratio - expected) / expected
        worst = max(worst, dev)
        msgs.append(f"  {name}: S_A/S_B={ratio:.6f} exp(-TE/T2)={expected:.6f} "
                    f"(TE={r.te*1e3:.1f}ms T2={t2*1e3:.1f}ms) dev={dev:.4f}")
    return worst <= 0.02, worst, msgs


# --- pytest entry points (validation-analytic tier; needs pytest) ------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    return binary


def test_t2_decay_cpu():
    """CPU: intra/extra b0 magnitude carries exp(-TE/T2) for each compartment."""
    ok, worst, msgs = evaluate(run_t2(_binary_or_skip(), n=6000, gpu=0))
    assert ok, f"T2 weighting off (worst rel dev {worst:.4f}):\n" + "\n".join(msgs)


def test_t2_decay_gpu():
    """GPU: same T2 weighting (self-skips without CUDA)."""
    import pytest
    r = run_t2(_binary_or_skip(), n=6000, gpu=1)
    if r.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    ok, worst, msgs = evaluate(r)
    assert ok, f"T2 weighting off (GPU, worst rel dev {worst:.4f}):\n" + "\n".join(msgs)


# --- T2 levels: exp(-TE/T2) at a low AND a high extra-compartment T2 ----------
T2_LEVELS = (0.040, 0.150)      # low / high (s), fixed diffusivity


def _assert_t2_levels(binary, gpu):
    import pytest
    for t2 in T2_LEVELS:
        r = run_t2(binary, n=6000, gpu=gpu, t2_extra_s=t2)
        if r.skipped:
            pytest.skip("binary built without CUDA / no GPU available")
        ok, worst, msgs = evaluate(r)
        assert ok, (f"T2={t2*1e3:.0f}ms weighting off "
                    f"({'GPU' if gpu else 'CPU'}, worst {worst:.4f}):\n" + "\n".join(msgs))


def test_t2_levels_cpu():
    """CPU: exp(-TE/T2) holds at both a low (40 ms) and a high (150 ms) T2."""
    _assert_t2_levels(_binary_or_skip(), gpu=0)


def test_t2_levels_gpu():
    """GPU: same (self-skips without CUDA)."""
    _assert_t2_levels(_binary_or_skip(), gpu=1)


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=6000)
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    ap.add_argument("--rtol", type=float, default=0.02,
                    help="max allowed relative deviation of the T2 ratio / fit")
    a = ap.parse_args()

    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1
    r = run_t2(binary, n=a.n, gpu=a.gpu)
    if r.skipped:
        print("SKIP: binary built without CUDA / no GPU available.")
        return 0

    ok, worst, msgs = evaluate(r)
    ok = worst <= a.rtol
    print(f"[T2 decay | {'multi-TE fit' if r.multi else 'single-TE differential'} | "
          f"N={a.n} {'GPU' if a.gpu else 'CPU'}]")
    for m in msgs:
        print(m)
    print(f"{'PASS' if ok else 'FAIL'}: worst rel dev={worst:.4f} (rtol {a.rtol})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
