#!/usr/bin/env python3
"""
Unit check for per-obstacle columns in spheres_list / cylinders_list.

Everything in a SI .conf (the default) is in standard units, so the T2 column of a sphere /
cylinder list is in SECONDS. An impermeable obstacle with T2 = 0.04 s, all walkers seeded
intra and a b0 acquisition at TE = 45 ms must give  S/N = exp(-TE/T2) = 0.325. (Read as
milliseconds -- the old behaviour -- the signal would be exp(-1125) ~ 0.)
"""
import math, os, subprocess, sys, tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "MC-DC_Simulator")
SCHEME = os.path.join(REPO, "tests", "accuracy", "PGSE_sample_scheme.scheme")  # row 0: b0, TE = 45 ms
N, TE, T2, TOL = 2000, 0.045, 0.04, 0.02


def run(kind, work):
    lst = os.path.join(work, kind + ".list")
    with open(lst, "w") as fh:
        fh.write("1e-6\n0 0 0 5 0 %g\n" % T2 if kind == "spheres" else "1e-6\n0 0 0 5 0 %g\n" % T2)
    prefix = os.path.join(work, kind)
    conf = prefix + ".conf"
    with open(conf, "w") as fh:
        fh.write("N %d\nT 200\nduration 0.045\ndiffusivity 0.6e-9\nt2_extra 1e9\nexp_prefix %s\n"
                 "scheme_file %s\nwrite_txt 1\nwrite_bin 0\nwrite_traj_file 0\ndeportation 0\n"
                 "ini_walkers_pos intra\n<obstacle>\n%s_list %s\n</obstacle>\n"
                 "<voxels>\n-1e-5 -1e-5 -1e-5\n1e-5 1e-5 1e-5\n</voxels>\nnum_process 1\nseed 7\n<END>\n"
                 % (N, prefix, SCHEME, kind[:-1] + "s" if kind == "spheres" else "cylinders", lst))
    r = subprocess.run([BIN, "--conf", conf], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0 or "[ERROR]" in r.stdout:
        sys.stderr.write(r.stdout); return None
    return float(open(prefix + "_DWI.txt").read().split()[0]) / N


def main():
    ok = True
    want = math.exp(-TE / T2)
    with tempfile.TemporaryDirectory() as work:
        for kind in ("spheres", "cylinders"):
            s = run(kind, work)
            good = s is not None and abs(s - want) <= TOL
            print("[%s_list] T2 column in seconds: S/N = %s, expected %.4f -> %s"
                  % (kind[:-1], "n/a" if s is None else "%.4f" % s, want, "PASS" if good else "FAIL"))
            ok &= good
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
