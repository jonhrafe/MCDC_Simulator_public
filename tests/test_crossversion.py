#!/usr/bin/env python3
"""
Cross-version check: the current build must reproduce the signal of an untouched historical release.

The old binary (v1.5 / v2.0) is built from this repo's own git history by tests/tools/build_old_version.py
(cached under tests/crossversion/old_binaries/, rebuilt only if missing). Both binaries run the same
experiment: the old one from its per-version config (tests/crossversion/confs/<version>/<case>.conf, old
unit convention), the current one from the SI config tests/accuracy/<case>.conf, CPU path, same N.

The two runs use different RNGs, so they are compared statistically. Signals are normalised by the mean of
the b=0 rows (E = S/S0). Per scheme row the Monte Carlo standard error of one run is sqrt((1 - E^2)/(2N))
(variance of cos(phase) for a Gaussian phase), and of the difference SE = sqrt(SE_old^2 + SE_cur^2).
FAIL if any row has |z| > --zmax (default 4.5: ~0.1 % false-alarm over 123 rows) or if a b-shell mean
(directions averaged) differs by more than --zshell standard errors (default 4).

Usage: test_crossversion.py --bin MC-DC_Simulator --version v2.0 --case cylinder_restricted [--n 20000]
Exit 0 = pass.
"""
import argparse, os, re, shutil, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "tools"))
import build_old_version  # noqa: E402

SCHEME = os.path.join(REPO, "tests", "accuracy", "PGSE_sample_scheme.scheme")


def rewrite_conf(src, dst, prefix, n):
    """copy a config, forcing exp_prefix, N, the CPU path and absolute repo paths."""
    out = []
    for line in open(src):
        key = line.split()[0].lower() if line.split() else ""
        if key == "exp_prefix":
            line = "exp_prefix %s\n" % prefix
        elif key == "n" and n:
            line = "N %d\n" % n
        elif key == "gpu":
            line = "gpu 0\n"
        else:
            line = re.sub(r"(?<![\w/])tests/", os.path.join(REPO, "tests") + "/", line)
        out.append(line)
    open(dst, "w").writelines(out)


def run(binary, conf, prefix, cwd):
    r = subprocess.run([binary, "--conf", conf], cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    sig = prefix + "_DWI.txt"
    if r.returncode != 0 or not os.path.exists(sig):
        sys.stdout.write(r.stdout[-3000:])
        raise RuntimeError("%s failed (rc=%d)" % (binary, r.returncode))
    return np.loadtxt(sig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=os.path.join(REPO, "MC-DC_Simulator"))
    ap.add_argument("--version", required=True, choices=("v1.5", "v2.0"))
    ap.add_argument("--case", required=True)
    ap.add_argument("--n", type=int, default=None, help="override N in both configs")
    ap.add_argument("--zmax", type=float, default=4.5)
    ap.add_argument("--zshell", type=float, default=4.0)
    a = ap.parse_args()

    old_conf = os.path.join(HERE, "crossversion", "confs", a.version, a.case + ".conf")
    cur_conf = os.path.join(HERE, "accuracy", a.case + ".conf")
    for c in (old_conf, cur_conf):
        if not os.path.exists(c):
            print("missing config: " + c); return 1
    old_bin = build_old_version.ensure_binary(a.version)

    work = tempfile.mkdtemp(prefix="mcdc_xv_")
    try:
        sigs, ns = {}, {}
        for tag, binary, conf in (("old", old_bin, old_conf), ("cur", os.path.abspath(a.bin), cur_conf)):
            prefix = os.path.join(work, tag)
            c = os.path.join(work, tag + ".conf")
            rewrite_conf(conf, c, prefix, a.n)
            ns[tag] = int(next(l.split()[1] for l in open(c) if l.split() and l.split()[0].lower() == "n"))
            sigs[tag] = run(str(binary), c, prefix, work)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    sch = np.loadtxt(SCHEME, skiprows=1)
    G = sch[:, 3]; b0 = G == 0
    E = {k: v / v[b0].mean() for k, v in sigs.items()}
    var = {k: np.clip(1 - E[k] ** 2, 1e-6, None) / 2 for k in E}
    se = np.sqrt(var["old"] / ns["old"] + var["cur"] / ns["cur"])
    d = E["cur"] - E["old"]; z = d / se
    ok = np.abs(z[~b0]).max() <= a.zmax

    # b-shells: unique (G, Delta, delta) among the diffusion-weighted rows
    key = np.round(sch[:, 3:6], 6)
    shells = sorted({tuple(k) for k, w in zip(key, ~b0) if w})
    print("%s vs current, case %s, N old=%d cur=%d" % (a.version, a.case, ns["old"], ns["cur"]))
    print("  rows: max|z| = %.2f (limit %.1f), rms diff = %.5f, mean z = %+.2f"
          % (np.abs(z[~b0]).max(), a.zmax, np.sqrt(np.mean(d[~b0] ** 2)), z[~b0].mean()))
    for s in shells:
        m = np.all(key == s, axis=1)
        dz = d[m].mean() / np.sqrt(np.sum(se[m] ** 2)) * m.sum()
        flag = abs(dz) > a.zshell
        ok &= not flag
        print("  shell G=%.4g Delta=%.4g delta=%.4g: E old %.4f cur %.4f, diff %+.5f, z %+.2f%s"
              % (s[0], s[1], s[2], E["old"][m].mean(), E["cur"][m].mean(), d[m].mean(), dz, "  <-- FAIL" if flag else ""))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
