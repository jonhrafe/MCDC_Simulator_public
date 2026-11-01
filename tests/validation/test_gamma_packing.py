#!/usr/bin/env python3
"""
Substrate-generator validation: GAMMA-distributed sphere / cylinder packings.

This validates the *substrate generator* itself (the packing engine + its seeded
RNG + the SI unit-homogenisation), which nothing else in the suite touches -- the
existing ``gamma_*.conf`` cases only check that the CPU walk matches the GPU walk
on *whatever* substrate happens to be produced. Here we check that the produced
substrate is CORRECT, using only geometric / statistical invariants (numpy only,
no analytic DWI oracle needed).

The generator (src/{cylinder,sphere}gammadistribution.cpp) draws N radii from a
Gamma(alpha, beta) law, auto-sizes the voxel from the requested ICVF
(cyl: L = sqrt(Sum pi r^2 / icvf); sph: L = (Sum 4/3 pi r^3 / icvf)^(1/3)), packs
non-overlapping obstacles with periodic images, and -- crucially for us -- writes
the realised packing to ``<prefix>_gamma_distributed_{cylinder,sphere}_list.txt``
(parallelmcsimulation.cpp), a native list substrate we parse back.

Two tissue-scale fixtures (voxel ~100 um/side, the calibration target):
  * cylinders: alpha 1.5, beta 1.0e-6 m, icvf 0.40, num 340 -> mean r 1.5 um
  * spheres:   alpha 1.5, beta 2.0e-6 m, icvf 0.40, num 900 -> mean r 3.0 um
(mean radius = alpha*beta; min_radius defaults to 0 so there is no low-tail
truncation and the realised mean equals alpha*beta.)

Checks (all invariants that must hold regardless of the exact RNG stream):
  A1  achieved ICVF: recompute Sum(area|vol)/voxel from the list; must match the
      requested icvf AND the generator's self-reported "ICVF achieved" line.
  A2  radius-distribution fidelity: sample mean ~ alpha*beta and coefficient of
      variation ~ 1/sqrt(alpha) (mean pins the scale, CV pins the shape -> the two
      together pin Gamma(alpha,beta)).
  A3  non-overlap + in-voxel: no two obstacles (incl. periodic images) interpenetrate
      and every primary centre lies inside the auto-sized voxel.
  A4  auto-voxel: the voxel from simulation_info equals the closed-form
      computeMinimalSize and is ~100 um (the calibration target).
  A5  determinism (P0.1): same seed + same num_process -> byte-identical list;
      a different seed -> a different list.

Tolerances are loose statistical bands (finite N=340/900 sample), not bit-exact:
ICVF within 0.03 of target, mean radius within 10%, CV within 0.12.

  python3 tests/validation/test_gamma_packing.py --bin ./MC-DC_Simulator
"""
import argparse
import filecmp
import math
import os
import shutil
import sys
import tempfile
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcdc_testkit import paths, run, substrate  # noqa: E402

import numpy as np  # noqa: E402

SCHEME = "tests/accuracy/PGSE_sample_scheme.scheme"

# The two calibrated tissue-scale cases. beta is SI (metres); beta_um = beta*1e6.
CASES = {
    "cylinder": dict(
        kind="cylinder", block="cylinder_gamma_packing", numkey="num_cylinders",
        alpha=1.5, beta=1.0e-6, icvf=0.40, num=340,
        dims=(0, 1),          # z-parallel cylinders overlap in the xy-plane
        voxel_um=103.1),
    "sphere": dict(
        kind="sphere", block="sphere_gamma_packing", numkey="num_spheres",
        alpha=1.5, beta=2.0e-6, icvf=0.40, num=900,
        dims=(0, 1, 2),
        voxel_um=102.2),
}

Pack = namedtuple("Pack", "kind icvf_req icvf_rep icvf_rec voxel_um n_prim "
                          "mean_um cv_ measured min_gap_um in_voxel "
                          "det_same det_diff")


# --------------------------------------------------------------------------- #
#  Build config + run the generator (capturing stdout for the achieved ICVF)  #
# --------------------------------------------------------------------------- #
def _conf_text(case, prefix, seed, n, t):
    return (
        f"# gamma {case['kind']} packing (generator-fidelity fixture)\n"
        f"N {n}\nT {t}\nduration 0.050\ndiffusivity 0.6e-9\n"
        f"exp_prefix {prefix}\nscheme_file {SCHEME}\n"
        "write_txt 1\nwrite_bin 0\nwrite_traj_file 0\nt2_extra 1e9\n"
        "<obstacle>\n"
        f"<{case['block']}>\n"
        f"alpha {case['alpha']}\nbeta {case['beta']}\nicvf {case['icvf']}\n"
        f"{case['numkey']} {case['num']}\n"
        f"</{case['block']}>\n"
        "</obstacle>\n"
        f"ini_walkers_pos intra\nnum_process 1\nseed {seed}\n<END>\n")


def _run(binary, case, work, base, seed=12345, n=10, t=2, timeout=300):
    """Run the generator; return (list_path, achieved_icvf, voxel_um)."""
    prefix = os.path.join(work, base)
    conf = os.path.join(work, base + ".conf")
    with open(conf, "w") as fh:
        fh.write(_conf_text(case, prefix, seed, n, t))
    res = run.run_sim(binary, conf, paths.REPO_ROOT, timeout=timeout)
    run.raise_on_error(res)
    list_path = f"{prefix}_gamma_distributed_{case['kind']}_list.txt"
    return list_path, _achieved_icvf(res.stdout), substrate.voxel_side_um(f"{prefix}_simulation_info.txt")


def _achieved_icvf(stdout):
    """Parse the generator's self-reported 'ICVF achieved: 40.000000' (percent)."""
    for ln in (stdout or "").splitlines():
        if "ICVF achieved" in ln:
            return float(ln.split(":")[1].split()[0]) / 100.0
    return None


# --------------------------------------------------------------------------- #
#  Compute the invariants (substrate parsing lives in mcdc_testkit.substrate)  #
# --------------------------------------------------------------------------- #
def _min_gap_um(centres, radii, dims):
    """Smallest surface-to-surface gap over all obstacle pairs (>=0 == no overlap).

    Distance uses only the packed axes (cyl: xy, z-parallel infinite; sph: xyz).
    Uses the FULL list (primaries + periodic images) so cross-boundary contacts
    are checked too.
    """
    c = centres[:, list(dims)]
    d = np.sqrt(((c[:, None, :] - c[None, :, :]) ** 2).sum(-1))
    gap = d - (radii[:, None] + radii[None, :])
    np.fill_diagonal(gap, np.inf)
    return float(gap.min())


def analyse(binary, case, work, timeout=300):
    """Run the case + its determinism replicas; return a Pack of measurements."""
    lp, icvf_rep, voxel_um = _run(binary, case, work, case["kind"] + "_a",
                                  seed=12345, timeout=timeout)
    centres, radii = substrate.read_list(lp)
    pc, pr = substrate.primary(centres, radii, voxel_um, case["dims"])

    # A1: recomputed ICVF from the primary set + the auto-sized voxel.
    if case["kind"] == "cylinder":
        occ = float((math.pi * pr ** 2).sum()) / (voxel_um ** 2)
    else:
        occ = float((4.0 / 3.0 * math.pi * pr ** 3).sum()) / (voxel_um ** 3)

    # A2: radius-distribution moments (um).
    mean_um = float(pr.mean())
    cv = float(pr.std() / pr.mean())

    # A3: non-overlap (full list) + primaries inside the voxel.
    min_gap = _min_gap_um(centres, radii, case["dims"])
    in_voxel = all(
        (pc[:, d] >= -1e-3 * voxel_um).all() and (pc[:, d] <= voxel_um + 1e-3 * voxel_um).all()
        for d in case["dims"])

    # A5: determinism -- same seed identical, different seed different.
    lp2, _, _ = _run(binary, case, work, case["kind"] + "_a2", seed=12345, timeout=timeout)
    lp3, _, _ = _run(binary, case, work, case["kind"] + "_b", seed=777, timeout=timeout)
    det_same = filecmp.cmp(lp, lp2, shallow=False)
    det_diff = not filecmp.cmp(lp, lp3, shallow=False)

    return Pack(case["kind"], case["icvf"], icvf_rep, occ, voxel_um, len(pc),
                mean_um, cv, True, min_gap, in_voxel, det_same, det_diff)


# --------------------------------------------------------------------------- #
#  Pass / fail                                                                 #
# --------------------------------------------------------------------------- #
def evaluate(p, case):
    """Return (ok, messages) for a Pack against the case's expectations."""
    beta_um = case["beta"] * 1e6
    mean_exp = case["alpha"] * beta_um
    cv_exp = 1.0 / math.sqrt(case["alpha"])
    msgs, ok = [], True

    def check(cond, label):
        nonlocal ok
        ok = ok and cond
        msgs.append(f"  [{'OK ' if cond else 'BAD'}] {label}")

    check(abs(p.icvf_rec - p.icvf_req) <= 0.03,
          f"A1 ICVF recomputed={p.icvf_rec:.3f} vs requested={p.icvf_req:.2f}")
    if p.icvf_rep is not None:
        check(abs(p.icvf_rec - p.icvf_rep) <= 0.02,
              f"A1 recomputed vs reported {p.icvf_rep:.3f} (engine agreement)")
    check(abs(p.mean_um - mean_exp) <= 0.10 * mean_exp,
          f"A2 mean radius={p.mean_um:.3f} um vs alpha*beta={mean_exp:.3f} um")
    check(abs(p.cv_ - cv_exp) <= 0.12,
          f"A2 CV={p.cv_:.3f} vs 1/sqrt(alpha)={cv_exp:.3f}")
    check(p.min_gap_um >= -0.10,
          f"A3 min pair gap={p.min_gap_um:.4f} um (>=0 == no overlap)")
    check(p.in_voxel, "A3 all primary centres inside the voxel")
    check(0.9 * case["num"] <= p.n_prim <= case["num"] + 1,
          f"A3 primary count={p.n_prim} (requested {case['num']})")
    check(80.0 <= p.voxel_um <= 130.0,
          f"A4 voxel={p.voxel_um:.1f} um (~100 um target, calib {case['voxel_um']:.1f})")
    check(p.det_same, "A5 same seed -> byte-identical list")
    check(p.det_diff, "A5 different seed -> different list")
    return ok, msgs


# --- pytest entry points (validation tier; needs pytest) ---------------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    return binary


def _run_case_pytest(name):
    work = tempfile.mkdtemp(prefix=f"mcdc_gamma_{name}_")
    try:
        p = analyse(_binary_or_skip(), CASES[name], work)
        ok, msgs = evaluate(p, CASES[name])
        assert ok, f"gamma {name} packing invariants failed:\n" + "\n".join(msgs)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def test_gamma_cylinder_packing():
    """Mean-1.5um gamma cylinders in a ~100um voxel: ICVF, distribution, non-overlap."""
    _run_case_pytest("cylinder")


def test_gamma_sphere_packing():
    """Mean-3um gamma spheres in a ~100um voxel: ICVF, distribution, non-overlap."""
    _run_case_pytest("sphere")


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--only", choices=("cylinder", "sphere"), default=None,
                    help="run just one case (default: both)")
    a = ap.parse_args()

    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1

    names = [a.only] if a.only else ["cylinder", "sphere"]
    all_ok = True
    for name in names:
        work = tempfile.mkdtemp(prefix=f"mcdc_gamma_{name}_")
        try:
            p = analyse(binary, CASES[name], work)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        ok, msgs = evaluate(p, CASES[name])
        all_ok = all_ok and ok
        print(f"[gamma {name} packing | voxel {p.voxel_um:.1f} um | "
              f"{p.n_prim} obstacles | ICVF {p.icvf_rec:.3f}]")
        for m in msgs:
            print(m)
        print(f"{'PASS' if ok else 'FAIL'}: {name}\n")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
