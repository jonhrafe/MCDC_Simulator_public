#!/usr/bin/env python3
"""
Analytic validation: restricted diffusion in a GAMMA PACKING vs the volume-weighted
Gaussian-Phase-Approximation (experiments 1 & 2).

Runs the calibrated tissue-scale packings (``gamma_{sphere,cylinder}_100um.conf``,
~100 µm voxel) with walkers seeded UNIFORMLY in the intra volume
(``ini_walkers_pos intra``) and the obstacles IMPERMEABLE, then compares the
Monte-Carlo intra signal to a closed-form oracle:

    S_ana(b) = Σ_i V_i · S_GPA(b; r_i) / Σ_i V_i

the **volume-weighted** single-obstacle GPA over the packing's own radii
(``mcdc_testkit.oracles.packing_signal``): spheres by V∝r³ (Neuman / Murday-Cotts,
isotropic), cylinders by V∝r² (van Gelderen, per gradient direction; all cylinders
lie along z). This is exact in expectation because uniform intra seeding puts a
walker count in each obstacle ∝ its volume — so the MC signal IS the volume-weighted
average of the per-obstacle restricted signals.

The GPA is the SAME model as dmipy's ``C4CylinderGaussianPhaseApproximation`` /
``S4SphereGaussianPhaseApproximation`` (see ``script_comparing_3v.py`` and
``tests/tools/gen_analytic.py``), re-implemented in pure numpy so the oracle is
evaluated live from each packing's radii — no dmipy, no frozen curve (the curve is
specific to the exact packing).

VALIDITY. GPA = infinite (cyl) / closed (sph), IMPERMEABLE, small b (≲5000 s/mm²).
The sample scheme tops out at b=2000 s/mm² → well inside the valid regime. Measured
agreement on this build: worst per-shell |MC − GPA| ≈ 0.0015 (sphere) / 0.0008 (cyl)
at N=30000 — the tolerance below (0.02) is a comfortable statistical band.

  python3 tests/validation/test_packing_restricted.py --bin ./MC-DC_Simulator --n 30000
  python3 tests/validation/test_packing_restricted.py --only sphere --plot report/
"""
import argparse
import os
import shutil
import sys
import tempfile
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcdc_testkit import io, metrics, oracles, paths, run, substrate  # noqa: E402
from mcdc_testkit import scheme as kscheme  # noqa: E402

import numpy as np  # noqa: E402

CASES = {
    "sphere":   dict(conf="gamma_sphere_100um.conf",   geometry="sphere",   dims=(0, 1, 2)),
    "cylinder": dict(conf="gamma_cylinder_100um.conf", geometry="cylinder", dims=(0, 1)),
}

Result = namedtuple("Result", "kind sc n_obs mc ana bvals skipped")


def run_case(binary, name, n=30000, gpu=0, timeout=1200):
    """Run the packing intra + compute the volume-weighted GPA over its radii."""
    case = CASES[name]
    conf_text = open(paths.accuracy(case["conf"])).read()
    D = float(io.read_conf_value(conf_text, "diffusivity"))
    sc = kscheme.parse_scheme(os.path.join(paths.REPO_ROOT,
                                           io.read_conf_value(conf_text, "scheme_file")))
    work = tempfile.mkdtemp(prefix=f"mcdc_pack_{name}_")
    try:
        prefix = run.run_case(binary, conf_text, work, name, paths.REPO_ROOT,
                              overrides={"N": n}, gpu=gpu, timeout=timeout)
        if prefix is None:
            return Result(name, sc, 0, None, None, sc.bvals, True)
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=n)
        mc = metrics.normalize_by_b0(dwi, sc.b0_mask)
        # radii the walkers actually saw: primary (de-duped) obstacles of THIS run
        voxel_um = substrate.voxel_side_um(prefix + "_simulation_info.txt")
        centres, radii = substrate.read_list(
            f"{prefix}_gamma_distributed_{name}_list.txt")
        _, pr = substrate.primary(centres, radii, voxel_um, case["dims"])
        radii_m = pr * 1e-6
        ana = oracles.packing_signal(sc, radii_m, D, case["geometry"])
        return Result(name, sc, len(pr), mc, ana, sc.bvals, False)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def evaluate(r, tol=0.02):
    """Per-shell |MC − GPA| check on the non-b0 shells. Returns (ok, worst, msgs)."""
    nb = ~r.sc.b0_mask
    mmc = metrics.per_shell_mean(r.mc[nb], r.bvals[nb])
    man = metrics.per_shell_mean(r.ana[nb], r.bvals[nb])
    worst, msgs = 0.0, []
    for b in sorted(mmc):
        d = abs(mmc[b] - man[b])
        worst = max(worst, d)
        msgs.append(f"  b={b:5.0f} s/mm²  MC={mmc[b]:.4f}  GPA_vw={man[b]:.4f}  |Δ|={d:.4f}")
    nrmse = metrics.nrmse(r.mc[nb], r.ana[nb])
    msgs.append(f"  NRMSE(full, per-measurement) = {nrmse:.4f}")
    return worst <= tol, worst, msgs


def write_plot(results, out_dir):
    """Per-shell overlay (per-shell scatter + mean-per-shell), matplotlib optional.

    Always writes a markdown data table; adds PNGs when matplotlib is importable.
    """
    os.makedirs(out_dir, exist_ok=True)
    md = ["# Packing restricted diffusion — MC vs volume-weighted GPA\n"]
    for r in results:
        if r.skipped:
            continue
        nb = ~r.sc.b0_mask
        mmc = metrics.per_shell_mean(r.mc[nb], r.bvals[nb])
        man = metrics.per_shell_mean(r.ana[nb], r.bvals[nb])
        md.append(f"\n## {r.kind} ({r.n_obs} obstacles)\n")
        md.append("| b (s/mm²) | MC | GPA (vol-weighted) | |Δ| |")
        md.append("|---|---|---|---|")
        for b in sorted(mmc):
            md.append(f"| {b:.0f} | {mmc[b]:.4f} | {man[b]:.4f} | {abs(mmc[b]-man[b]):.4f} |")
    with open(os.path.join(out_dir, "packing_restricted.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    for r in results:
        if r.skipped:
            continue
        sc = r.sc
        nb = ~sc.b0_mask
        ub = metrics.shells(sc.bvals[nb])
        xang = 90 - kscheme.angle_to_axis_deg(sc.gdirs)      # 0=perp, 90=|| z
        fig, axs = plt.subplots(1, len(ub), figsize=(5 * len(ub), 4), squeeze=False)
        for j, b in enumerate(ub):
            m = nb & (np.round(sc.bvals / 1e6, 0) == b)
            ax = axs[0][j]
            o = np.argsort(xang[m])
            ax.scatter(xang[m][o], r.mc[m][o], marker="x", color="tab:red", label="MC")
            ax.plot(xang[m][o], r.ana[m][o], "-k", lw=2, label="GPA (vw)")
            ax.set_title(f"{r.kind}  b={b:.0f}")
            ax.set_xlabel("90 − angle to z [deg]")
            ax.set_ylabel("S/S0")
            ax.grid(True, ls=":", alpha=0.4)
            ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"packing_{r.kind}.png"), dpi=120)
        plt.close(fig)


# --- pytest entry points (analytic tier; needs pytest) -----------------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    return binary


def _vs_gpa(name, gpu):
    """Run the packing (CPU or GPU) and check vs the GPA; skip if the GPU is absent."""
    import pytest
    r = run_case(_binary_or_skip(), name, n=30000, gpu=gpu)
    if r.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    ok, worst, msgs = evaluate(r)
    assert ok, f"{name} packing ({'GPU' if gpu else 'CPU'}) off GPA " \
               f"(worst {worst:.4f}):\n" + "\n".join(msgs)


def test_packing_sphere_restricted():
    """CPU: sphere packing intra ≈ volume-weighted Neuman GPA (experiment 1)."""
    _vs_gpa("sphere", gpu=0)


def test_packing_cylinder_restricted():
    """CPU: cylinder packing intra ≈ area-weighted van Gelderen GPA (experiment 2)."""
    _vs_gpa("cylinder", gpu=0)


def test_packing_sphere_restricted_gpu():
    """GPU: same as experiment 1 on the CUDA path (self-skips without CUDA)."""
    _vs_gpa("sphere", gpu=1)


def test_packing_cylinder_restricted_gpu():
    """GPU: same as experiment 2 on the CUDA path (self-skips without CUDA)."""
    _vs_gpa("cylinder", gpu=1)


def _parity(name):
    """CPU ≡ GPU on the identical packing: per-shell means agree within MC noise."""
    import pytest
    binary = _binary_or_skip()
    g = run_case(binary, name, n=30000, gpu=1)
    if g.skipped:
        pytest.skip("binary built without CUDA / no GPU available")
    c = run_case(binary, name, n=30000, gpu=0)
    worst = metrics.per_shell_max_dev(c.mc[~c.sc.b0_mask], g.mc[~g.sc.b0_mask],
                                      c.bvals[~c.sc.b0_mask])
    assert worst <= 0.02, f"{name} packing CPU/GPU parity off: worst per-shell |Δ|={worst:.4f}"


def test_packing_sphere_cpu_gpu_parity():
    """Experiment 1: CPU ≡ GPU (same answer)."""
    _parity("sphere")


def test_packing_cylinder_cpu_gpu_parity():
    """Experiment 2: CPU ≡ GPU (same answer)."""
    _parity("cylinder")


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=30000)
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--only", choices=("sphere", "cylinder"), default=None)
    ap.add_argument("--plot", default=None, metavar="DIR",
                    help="write a markdown table (+ PNGs if matplotlib) to DIR")
    a = ap.parse_args()

    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1

    names = [a.only] if a.only else ["sphere", "cylinder"]
    results, all_ok = [], True
    for name in names:
        r = run_case(binary, name, n=a.n, gpu=a.gpu)
        results.append(r)
        if r.skipped:
            print(f"SKIP {name}: binary built without CUDA / no GPU.")
            continue
        ok, worst, msgs = evaluate(r, tol=a.tol)
        all_ok = all_ok and ok
        print(f"[{name} packing intra | {r.n_obs} obstacles | N={a.n} "
              f"{'GPU' if a.gpu else 'CPU'}]")
        for m in msgs:
            print(m)
        print(f"{'PASS' if ok else 'FAIL'}: worst |MC−GPA|={worst:.4f} (tol {a.tol})\n")
    if a.plot:
        write_plot(results, a.plot)
        print(f"wrote report to {a.plot}/")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
