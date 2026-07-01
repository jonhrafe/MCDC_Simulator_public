#!/usr/bin/env python3
"""
Per-obstacle (list-column) permeability regression.

A spheres_list line is "x y z r kappa T2": the 5th column is a PER-SPHERE
permeability. This must actually drive membrane exchange. It regressed silently
because Sphere/Cylinder copy constructors did not copy the Obstacle base, so the
per-list kappa/T2/d_intra were reset to defaults on push_back (a per-sphere kappa
column then behaved fully impermeable; only the GLOBAL `permeability` worked,
because it is re-applied after insertion).

This test runs three otherwise-identical CPU simulations of one sphere:
  A) global `permeability`  + list kappa column = 0   (the working path)
  B) NO global permeability + list kappa column = K    (the per-obstacle path)
  C) NO global permeability + list kappa column = 0    (impermeable control)
and asserts:
  * B differs from C  -> the per-sphere column produces exchange (not dropped);
  * B matches  A      -> the per-sphere column gives the SAME exchange as global.

If the copy-ctor drop ever returns, B collapses onto C and this test fails.
CPU-only (gpu 0); self-contained (writes its own configs + lists to a temp dir).
"""
import glob
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KAPPA = 5e-4            # m/s; strong enough to give a clear exchange signature
# B-vs-C must exceed this (exchange present); B-vs-A must stay under it (matches global).
# The measured exchange here is ~16 pp of S/S0; MC noise at N is ~1-2 pp.
TOL = 0.05


def write_scheme(path):
    # Self-contained STEJSKALTANNER PGSE scheme (Gx Gy Gz |G| Delta delta TE): a b0 plus a
    # z-gradient ramp up to b ~ 2000 s/mm^2, where restricted-vs-permeable contrast is clear.
    lines = ["VERSION: STEJSKALTANNER", "0 0 0 0 0.025 0.005 0.045"]
    for G in (0.03, 0.06, 0.09, 0.12, 0.15, 0.155):
        lines.append("0 0 1 %f 0.025 0.005 0.045" % G)
    open(path, "w").write("\n".join(lines) + "\n")


def write_list(path, kappa):
    # scale 1e-6 m/unit; one sphere r=4 um at the origin, per-sphere kappa, T2 ~ inf.
    open(path, "w").write("1e-06\n0 0 0 4 %g 1000000000\n" % kappa)


def write_conf(path, prefix, list_file, global_perm, scheme):
    perm = ("permeability %g\n" % global_perm) if global_perm > 0 else ""
    open(path, "w").write(
        "N 10000\nT 1000\nduration 0.050\ndiffusivity 0.6e-9\n"
        "%sexp_prefix %s\nscheme_file %s\n"
        "write_txt 1\nwrite_bin 0\nt2_extra 1e9\nseparate_signals 1\ndeportation 0\n"
        "<obstacle>\nspheres_list %s\n</obstacle>\n"
        "<voxels>\n-5e-6 -5e-6 -5e-6\n 5e-6  5e-6  5e-6\n</voxels>\n"
        "num_process 1\nseed 12345\ngpu 0\n<END>\n"
        % (perm, prefix, scheme, list_file))


def run(binp, conf):
    res = subprocess.run([binp, "--conf", conf], cwd=REPO,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    if res.returncode != 0 or "[ERROR]" in (res.stdout or ""):
        sys.stderr.write(res.stdout or ""); return None
    prefix = conf[:-5]
    cands = [f for f in glob.glob(prefix + "*_DWI.txt") if not f.endswith(("_intra.txt", "_extra.txt"))]
    if not cands:
        return None
    f = max(cands, key=os.path.getmtime)
    return [float(x) for x in open(f) if x.strip()]


def norm_maxdiff(x, y):
    n = min(len(x), len(y)); s0 = max(max(x[:n]), max(y[:n]))
    return max(abs(x[i] - y[i]) for i in range(n)) / s0


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    work = tempfile.mkdtemp(prefix="mcdc_permlist_")
    try:
        scheme = os.path.join(work, "pgse.scheme"); write_scheme(scheme)
        write_list(os.path.join(work, "kappa.list"), KAPPA)
        write_list(os.path.join(work, "zero.list"), 0.0)
        write_conf(os.path.join(work, "A.conf"), os.path.join(work, "A"), os.path.join(work, "zero.list"),  KAPPA, scheme)  # global
        write_conf(os.path.join(work, "B.conf"), os.path.join(work, "B"), os.path.join(work, "kappa.list"), 0.0,   scheme)  # per-sphere
        write_conf(os.path.join(work, "C.conf"), os.path.join(work, "C"), os.path.join(work, "zero.list"),  0.0,   scheme)  # impermeable
        A = run(binp, os.path.join(work, "A.conf"))
        B = run(binp, os.path.join(work, "B.conf"))
        C = run(binp, os.path.join(work, "C.conf"))
        if A is None or B is None or C is None:
            print("FAIL: a simulation did not produce a DWI output."); return 1
        d_BC = norm_maxdiff(B, C)   # per-sphere vs impermeable -> must be LARGE (exchange present)
        d_BA = norm_maxdiff(B, A)   # per-sphere vs global      -> must be SMALL (same exchange)
        print("[permeability_list] per-sphere-vs-impermeable |dS/S0|=%.4f (need > %.2f), "
              "per-sphere-vs-global |dS/S0|=%.4f (need < %.2f)" % (d_BC, TOL, d_BA, TOL))
        if d_BC <= TOL:
            print("FAIL: per-sphere kappa column produced no exchange (copy-ctor drop regressed?)."); return 1
        if d_BA >= TOL:
            print("FAIL: per-sphere kappa disagrees with the same global kappa."); return 1
        print("PASS: per-sphere permeability column drives exchange and matches global kappa.")
        return 0
    finally:
        import shutil; shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
