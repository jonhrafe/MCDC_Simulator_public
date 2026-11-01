#!/usr/bin/env python3
"""
Time-dependence: restricted diffusion vs Δ, and permeability × Δ (the NEXI signature).

Uses ``tests/accuracy/nexi_time_dependence.scheme`` — a subset of an 8-Δ NEXI protocol:
Δ ∈ {10,20,30,45} ms, b ∈ {1000,2000,3000,5000,10000} s/mm², δ=4 ms, TE=58 ms (fixed).
CRITICAL: the conf must set ``duration ≥ 0.060`` (TE=58 ms; the engine has NO duration≥TE
guard, so a too-short duration silently truncates the second gradient lobe).

Two things, on one 5 µm sphere (native list AND a PLY mesh — same geometry, so one oracle):

  (a) RESTRICTION × Δ (impermeable): the signal matches Neuman GPA at each Δ for b ≤ 3000
      (GPA is valid only at small b), and the ADC DECREASES with Δ — restriction becomes
      visible as the diffusion time grows.

  (b) PERMEABILITY × Δ: with a permeable membrane, the exchange effect
      E(κ,Δ) = S_impermeable − S_permeable (same seed) grows with κ AND with Δ — permeability
      is nearly invisible at short Δ and clearly visible at long Δ (the NEXI idea). κ=0
      recovers (a). The crossed fraction f_extra grows with κ.

VALIDITY: the GPA is impermeable, small-b (≲5000). GPA-match runs only at b ≤ 3000; the
b=5000/10000 shells drive the CPU/GPU parity and the exchange-monotonicity checks (which
need no closed form). The quantitative per-encounter crossing probability is validated
separately in tests/test_permeability_powles.py; here the focus is the Δ-dependence.

  python3 tests/validation/test_time_dependence.py --bin ./MC-DC_Simulator
"""
import argparse
import os
import shutil
import sys
import tempfile
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcdc_testkit import io, metrics, oracles, paths, run  # noqa: E402
from mcdc_testkit import scheme as kscheme  # noqa: E402

import numpy as np  # noqa: E402

SCHEME = "tests/accuracy/nexi_time_dependence.scheme"
R = 5e-6                     # sphere radius (both native list and the scaled PLY)
D = 0.6e-9                   # intra diffusivity (m^2/s)
D_EXTRA = 2.0e-9             # extra diffusivity
DURATION = 0.060            # >= TE (0.058) -- REQUIRED, no engine guard
KAPPA = [0.0, 1e-5, 5e-5, 2e-4]      # none / low / med / high (m/s)
GPA_BMAX = 3000.0           # s/mm^2 -- GPA validity ceiling for the analytic match
NATIVE_LIST = "tests/accuracy/single_sphere.list"       # 5 um, kappa=0, T2 inf
PLY = "tests/accuracy/meshes/unitMesh.ply"              # ~1 um sphere at scale 0.001
PLY_SCALE = 0.005                                       # -> 5 um sphere (matches native)

Run = namedtuple("Run", "ss0 f_extra skipped")


def _abspath(rel):
    return os.path.join(paths.REPO_ROOT, rel)


def _run(binary, kind, kappa, n=15000, t=1000, gpu=0, timeout=400):
    """Run the 5 um sphere (kind='native'|'ply') at permeability kappa; return a Run.

    Both substrates are the SAME 5 um / D sphere: native via single_sphere.list (+ global
    obstacle_permeability), PLY via a per-kappa ply_extended_file_list on unitMesh scaled x5.
    """
    sc = kscheme.parse_scheme(_abspath(SCHEME))
    work = tempfile.mkdtemp(prefix=f"mcdc_td_{kind}_")
    try:
        perm_line = ""
        if kind == "ply":
            lst = os.path.join(work, "sphere.list")
            # "<ply> <scale> <d_intra> <t2 (s, must be >0)> <permeability>"; T2 huge -> no decay
            with open(lst, "w") as fh:
                fh.write(f"{_abspath(PLY)} {PLY_SCALE} {D} 1e9 {kappa}\n")
            obstacle = f"ply_extended_file_list {lst}"
        else:
            obstacle = f"spheres_list {_abspath(NATIVE_LIST)}"
            if kappa > 0:
                perm_line = f"obstacle_permeability {kappa}\n"   # TOP-LEVEL, not in <obstacle>
        prefix = os.path.join(work, "s")
        conf = (
            f"N {n}\nT {t}\nduration {DURATION}\ndiffusivity {D}\nd_extra {D_EXTRA}\n"
            f"{perm_line}"
            f"exp_prefix {prefix}\nscheme_file {SCHEME}\n"
            "write_txt 1\nwrite_bin 0\nwrite_traj_file 0\nseparate_signals\nt2_extra 1e9\n"
            f"<obstacle>\n{obstacle}\n</obstacle>\n"
            f"<voxels>\n{-1.2*R} {-1.2*R} {-1.2*R}\n{1.2*R} {1.2*R} {1.2*R}\n</voxels>\n"
            f"ini_walkers_pos intra\nnum_process 1\nseed 12345\ngpu {int(gpu)}\n<END>\n")
        conf_path = os.path.join(work, "c.conf")
        with open(conf_path, "w") as fh:
            fh.write(conf)
        res = run.run_sim(binary, conf_path, paths.REPO_ROOT, timeout=timeout)
        if gpu and run.is_gpu_skip(res.stdout or ""):
            return Run(None, None, True)
        run.raise_on_error(res)
        tot = io.load_dwi(prefix + "_DWI.txt", n_walkers=n)
        extra = io.load_dwi(prefix + "_DWI_extra.txt", n_walkers=n)
        b0 = sc.b0_mask
        f_extra = float(extra[b0].mean() / tot[b0].mean())
        return Run(metrics.normalize_by_b0(tot, b0), f_extra, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _by(sc, delta_ms, bmax=None):
    """Boolean mask for a given Delta (ms) and optional b<=bmax (s/mm^2), non-b0."""
    m = (~sc.b0_mask) & (np.round(sc.Delta * 1e3).astype(int) == delta_ms)
    if bmax is not None:
        m &= (sc.bvals / 1e6 <= bmax + 1e-6)
    return m


# --------------------------------------------------------------------------- #
#  Checks                                                                      #
# --------------------------------------------------------------------------- #
DELTAS = [10, 20, 30, 45]


def check_restriction(imp, sc, tol=0.03):
    """(a) impermeable sphere == Neuman GPA(Δ) at b<=3000, and ADC decreases with Δ."""
    gpa = oracles.sphere_neuman(sc, R, D)
    msgs, ok = [], True
    worst = 0.0
    for dl in DELTAS:
        m = _by(sc, dl, bmax=GPA_BMAX)
        w = float(np.max(np.abs(imp.ss0[m] - gpa[m])))
        worst = max(worst, w)
    ok = ok and worst <= tol
    msgs.append(f"  [{'OK ' if worst <= tol else 'BAD'}] κ=0 vs Neuman GPA (b≤3000): worst |Δ|={worst:.4f} (tol {tol})")
    # ADC at the lowest b vs Delta -> must decrease (restriction grows with time)
    adcs = []
    for dl in DELTAS:
        m = _by(sc, dl) & (np.round(sc.bvals / 1e6).astype(int) == 1000)
        adcs.append(float(-np.log(max(imp.ss0[m][0], 1e-9)) / (1000 * 1e6)))
    mono = all(a > b - 1e-13 for a, b in zip(adcs, adcs[1:]))
    ok = ok and mono
    msgs.append(f"  [{'OK ' if mono else 'BAD'}] ADC ↓ with Δ: "
                + ", ".join(f"{d}ms:{a*1e9:.3f}" for d, a in zip(DELTAS, adcs)) + " (1e-9)")
    return ok, msgs


def check_permeability(runs, sc, tol=0.03):
    """(b) exchange effect E(κ,Δ) grows with κ and Δ; f_extra grows with κ. runs: {κ: Run}."""
    imp = runs[0.0]
    msgs, ok = [], True

    # κ=0 must still match GPA (baseline integrity)
    ok0, _ = check_restriction(imp, sc, tol=tol)

    # E(κ,Δ) = mean_{b<=3000}(S_imp - S_perm)
    def E(kappa, dl):
        m = _by(sc, dl, bmax=GPA_BMAX)
        return float((imp.ss0[m] - runs[kappa].ss0[m]).mean())

    fes = [runs[k].f_extra for k in KAPPA]
    fe_mono = all(a < b + 1e-6 for a, b in zip(fes, fes[1:]))
    ok = ok and fe_mono
    msgs.append(f"  [{'OK ' if fe_mono else 'BAD'}] f_extra ↑ with κ: "
                + ", ".join(f"{k:g}:{f:.3f}" for k, f in zip(KAPPA, fes)))

    for kappa in KAPPA[1:]:
        Es = [E(kappa, dl) for dl in DELTAS]
        grows_delta = all(a < b + 1e-4 for a, b in zip(Es, Es[1:]))
        ok = ok and grows_delta
        msgs.append(f"  [{'OK ' if grows_delta else 'BAD'}] κ={kappa:g}: E(Δ) ↑ with Δ: "
                    + ", ".join(f"{d}ms:{e:.3f}" for d, e in zip(DELTAS, Es)))
    # E grows with κ at the longest Δ
    E_long = [E(k, DELTAS[-1]) for k in KAPPA[1:]]
    grows_kappa = all(a < b + 1e-4 for a, b in zip(E_long, E_long[1:]))
    ok = ok and grows_kappa
    msgs.append(f"  [{'OK ' if grows_kappa else 'BAD'}] E ↑ with κ (Δ=45): "
                + ", ".join(f"{k:g}:{e:.3f}" for k, e in zip(KAPPA[1:], E_long)))
    return ok and ok0, msgs


def _parity(cpu, gpu, sc):
    return metrics.per_shell_max_dev(cpu.ss0[~sc.b0_mask], gpu.ss0[~sc.b0_mask],
                                     sc.bvals[~sc.b0_mask])


# --- pytest entry points -----------------------------------------------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    if not os.path.exists(_abspath(SCHEME)):
        pytest.skip(f"scheme missing: {SCHEME}")
    return binary


def test_restriction_time_dependence():
    """(a) impermeable 5µm sphere: signal == Neuman GPA(Δ), ADC ↓ with Δ."""
    sc = kscheme.parse_scheme(_abspath(SCHEME))
    ok, msgs = check_restriction(_run(_binary_or_skip(), "native", 0.0, gpu=0), sc)
    assert ok, "restriction time-dependence failed:\n" + "\n".join(msgs)


def _perm_case(kind):
    sc = kscheme.parse_scheme(_abspath(SCHEME))
    binary = _binary_or_skip()
    runs = {k: _run(binary, kind, k, gpu=0) for k in KAPPA}
    ok, msgs = check_permeability(runs, sc)
    assert ok, f"permeability×Δ ({kind}) failed:\n" + "\n".join(msgs)


def test_permeability_time_dependence_native():
    """(b) native 5µm sphere: exchange effect grows with κ and Δ (NEXI signature)."""
    _perm_case("native")


def test_permeability_time_dependence_ply():
    """(b) PLY 5µm sphere: same, on the mesh permeability path."""
    _perm_case("ply")


def test_time_dependence_cpu_gpu_parity():
    """CPU ≡ GPU on the identical time-dependence run (native + PLY), self-skip w/o CUDA."""
    import pytest
    sc = kscheme.parse_scheme(_abspath(SCHEME))
    binary = _binary_or_skip()
    for kind in ("native", "ply"):
        g = _run(binary, kind, 5e-5, gpu=1)
        if g.skipped:
            pytest.skip("binary built without CUDA / no GPU available")
        c = _run(binary, kind, 5e-5, gpu=0)
        worst = _parity(c, g, sc)
        assert worst <= 0.02, f"{kind} CPU/GPU parity off: worst per-shell |Δ|={worst:.4f}"


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=15000)
    ap.add_argument("--kind", choices=("native", "ply", "both"), default="both")
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    a = ap.parse_args()
    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: binary not found: {binary}", file=sys.stderr)
        return 1
    sc = kscheme.parse_scheme(_abspath(SCHEME))
    kinds = ["native", "ply"] if a.kind == "both" else [a.kind]
    all_ok = True
    for kind in kinds:
        runs = {k: _run(binary, kind, k, n=a.n, gpu=a.gpu) for k in KAPPA}
        if runs[0.0].skipped:
            print(f"SKIP {kind}: no GPU."); continue
        rok, rmsgs = check_restriction(runs[0.0], sc)
        pok, pmsgs = check_permeability(runs, sc)
        all_ok = all_ok and rok and pok
        print(f"[time-dependence | {kind} sphere | N={a.n} {'GPU' if a.gpu else 'CPU'}]")
        for m in rmsgs + pmsgs:
            print(m)
        print(f"{'PASS' if rok and pok else 'FAIL'}: {kind}\n")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
