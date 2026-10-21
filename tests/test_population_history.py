#!/usr/bin/env python3
"""
Population-history and radii-injection tests (CPU, self-contained: writes its own configs/lists to a temp dir).

A) Four-population dynamical-history split (`crossed_filter`).
   One permeable sphere, uniform seeding, T2 = inf, num_process 1 + fixed seed (so the three
   runs below walk IDENTICAL trajectories; only the signal accumulation differs):
     crossed_filter 0  legacy: intra/extra split by FINAL location, all walkers
     crossed_filter 1  crossed-only,   split by SEED compartment (intra->extra / extra->intra)
     crossed_filter 2  uncrossed-only, split by SEED compartment (intra-never / hindered)
   Asserts
     * partition: total(cf0) == (intra+extra)(cf1) + (intra+extra)(cf2)  for every acquisition;
     * permeable: both efflux (intra seed) and influx (extra seed) channels are populated;
     * impermeable control: crossed channels are exactly empty, and the uncrossed channels
       reproduce the legacy total and the legacy intra signal.

B) Explicit radii injection (`radii_file` in <sphere_gamma_packing>).
   Asserts, from the emitted sphere list (x y z r perm T2 population_label is_replica):
     * the physical originals (is_replica == 0) are EXACTLY the injected (radius, label) multiset;
     * every periodic replica (is_replica == 1) copies an original's (radius, label);
     * (the packer keeps its best partial packing, so the target ICVF is set low enough that all
       injected spheres fit; a list that does not fit is truncated by the existing packer);
     * without radii_file, labels stay -1 and the original count equals num_spheres.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "MC-DC_Simulator")
REL_TOL = 1e-5          # partition identity (same trajectories; limited by the 6-digit text output)
LEAK_TOL = 2e-3         # impermeable control vs legacy: allows the ~0.1% membrane-leak walkers
KAPPA = 2e-5            # m/s: weak enough that all four channels stay populated on a 4 um sphere


def write_scheme(path):
    lines = ["VERSION: STEJSKALTANNER", "0 0 0 0 0.025 0.005 0.045"]
    for G in (0.03, 0.06, 0.09, 0.12, 0.15):
        lines.append("0 0 1 %f 0.025 0.005 0.045" % G)
    open(path, "w").write("\n".join(lines) + "\n")


def run(conf):
    res = subprocess.run([BIN, "--conf", conf], cwd=REPO, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, timeout=600)
    if res.returncode != 0 or "[ERROR]" in res.stdout:
        sys.stderr.write(res.stdout)
        raise RuntimeError("simulation failed: " + conf)
    return res.stdout


def latest(prefix, suffix):
    c = sorted(glob.glob(prefix + "*" + suffix), key=os.path.getmtime)
    if not c:
        raise RuntimeError("no output %s*%s" % (prefix, suffix))
    return c[-1]


def read_vec(path):
    return [float(x) for x in open(path) if x.strip()]


def signals(prefix):
    """(total, intra, extra) DWI vectors from the latest output of a run."""
    tot = read_vec(latest(prefix, "_DWI.txt"))
    intra = read_vec(latest(prefix, "_DWI_intra.txt"))
    extra = read_vec(latest(prefix, "_DWI_extra.txt"))
    return tot, intra, extra


# ----------------------------------------------------------------------------- A
def sphere_conf(path, prefix, scheme, lst, cf):
    flt = ("crossed_filter %d\n" % cf) if cf is not None else "separate_signals 1\n"
    open(path, "w").write(
        "N 8000\nT 1000\nduration 0.050\ndiffusivity 0.6e-9\n"
        "exp_prefix %s\nscheme_file %s\nwrite_txt 1\nwrite_bin 0\nt2_extra 1e9\n%s"
        "deportation 0\n<obstacle>\nspheres_list %s\n</obstacle>\n"
        "<voxels>\n-5e-6 -5e-6 -5e-6\n 5e-6  5e-6  5e-6\n</voxels>\n"
        "num_process 1\nseed 12345\n<END>\n" % (prefix, scheme, flt, lst))


def four_population(work, scheme):
    ok = True
    results = {}
    for tag, kappa in (("perm", KAPPA), ("imperm", 0.0)):
        lst = os.path.join(work, tag + ".list")
        open(lst, "w").write("1e-06\n0 0 0 4 %g 1000000000\n" % kappa)
        sig = {}
        for cf in (None, 1, 2):
            name = "%s_cf%s" % (tag, "0" if cf is None else cf)
            prefix = os.path.join(work, name)
            conf = prefix + ".conf"
            sphere_conf(conf, prefix, scheme, lst, cf)
            run(conf)
            sig[0 if cf is None else cf] = signals(prefix)
        results[tag] = sig

    # ---- partition identity (both substrates)
    for tag, sig in results.items():
        tot0 = sig[0][0]
        part = [sum(x) for x in zip(sig[1][1], sig[1][2], sig[2][1], sig[2][2])]
        err = max(abs(a - b) for a, b in zip(tot0, part)) / max(tot0)
        print("[four_pop/%s] partition |total - (cf1+cf2)|/S0 = %.2e (need < %.0e)" % (tag, err, REL_TOL))
        if err >= REL_TOL:
            print("FAIL: crossed/uncrossed channels do not partition the walkers (%s)." % tag)
            ok = False

    # ---- permeable: all four channels populated (b0 = walker counts, T2 = inf)
    p = results["perm"]
    ch = {"intra-never": p[2][1][0], "hindered": p[2][2][0],
          "efflux i->e": p[1][1][0], "influx e->i": p[1][2][0]}
    print("[four_pop/perm] b0 channel counts: " + ", ".join("%s=%.0f" % kv for kv in ch.items()))
    for k, v in ch.items():
        if v <= 0:
            print("FAIL: channel '%s' is empty on a permeable substrate." % k)
            ok = False

    # ---- impermeable control
    q = results["imperm"]
    crossed = max(abs(v) for v in q[1][1] + q[1][2])
    d_tot = max(abs(a - b) for a, b in zip(q[0][0], [x + y for x, y in zip(q[2][1], q[2][2])])) / max(q[0][0])
    d_intra = max(abs(a - b) for a, b in zip(q[0][1], q[2][1])) / max(q[0][0])
    print("[four_pop/imperm] max|crossed channel| = %.3g (need 0), uncrossed-vs-legacy total %.2e "
          "(need < %.0e), intra %.2e (need < %.0e)" % (crossed, d_tot, REL_TOL, d_intra, LEAK_TOL))
    if crossed != 0.0:
        print("FAIL: crossed channels non-empty without permeability."); ok = False
    if d_tot >= REL_TOL or d_intra >= LEAK_TOL:
        print("FAIL: uncrossed channels do not reproduce the legacy signal without permeability."); ok = False
    return ok


# ----------------------------------------------------------------------------- B
def gamma_conf(path, prefix, scheme, radii_file):
    rf = ("radii_file %s\n" % radii_file) if radii_file else ""
    open(path, "w").write(
        "N 50\nT 50\nduration 0.050\ndiffusivity 0.6e-9\n"
        "exp_prefix %s\nscheme_file %s\nwrite_txt 1\nwrite_bin 0\nt2_extra 1e9\n"
        "<obstacle>\n<sphere_gamma_packing>\nalpha 1.5\nbeta 0.5e-6\nicvf 0.08\nnum_spheres 12\n%s"
        "</sphere_gamma_packing>\n</obstacle>\n"
        "ini_walkers_pos intra\nnum_process 1\nseed 12345\n<END>\n" % (prefix, scheme, rf))


def read_sphere_list(prefix):
    f = latest(prefix, "_gamma_distributed_sphere_list.txt")
    rows = [l.split() for l in open(f).read().splitlines()[1:] if l.strip()]
    # x y z r perm T2 label is_replica   (r in um)
    return [(round(float(r[3]), 6), int(r[6]), int(r[7])) for r in rows]


def radii_injection(work, scheme):
    ok = True
    injected = [(1.5, 1), (1.5, 1), (1.5, 1), (1.0, 2), (1.0, 2), (1.0, 2),
                (0.8, 3), (0.8, 3), (0.8, 3), (0.6, 4), (0.6, 4), (0.6, 4)]
    rfile = os.path.join(work, "radii.txt")
    open(rfile, "w").write("".join("%g %d\n" % rl for rl in injected))
    prefix = os.path.join(work, "inj")
    gamma_conf(prefix + ".conf", prefix, scheme, rfile)
    run(prefix + ".conf")
    rows = read_sphere_list(prefix)
    orig = Counter((r, l) for r, l, rep in rows if rep == 0)
    repl = [(r, l) for r, l, rep in rows if rep == 1]
    want = Counter((round(r, 6), l) for r, l in injected)
    print("[radii_file] originals=%d replicas=%d (injected %d)" % (sum(orig.values()), len(repl), len(injected)))
    if orig != want:
        print("FAIL: originals %s != injected %s" % (dict(orig), dict(want))); ok = False
    bad = [rl for rl in repl if rl not in orig]
    if bad:
        print("FAIL: replicas that do not copy an original's (radius,label): %s" % bad[:5]); ok = False

    # control: no radii_file -> plain gamma, labels -1, original count == num_spheres
    prefix = os.path.join(work, "gam")
    gamma_conf(prefix + ".conf", prefix, scheme, None)
    run(prefix + ".conf")
    rows = read_sphere_list(prefix)
    n_orig = sum(1 for _, _, rep in rows if rep == 0)
    labels = set(l for _, l, _ in rows)
    print("[radii_file/control] originals=%d (num_spheres 12), labels=%s" % (n_orig, sorted(labels)))
    if n_orig != 12 or labels != {-1}:
        print("FAIL: gamma control changed (originals %d, labels %s)." % (n_orig, labels)); ok = False
    return ok


def main():
    work = tempfile.mkdtemp(prefix="mcdc_pophist_")
    try:
        scheme = os.path.join(work, "pgse.scheme"); write_scheme(scheme)
        ok_a = four_population(work, scheme)
        ok_b = radii_injection(work, scheme)
        if ok_a and ok_b:
            print("PASS: four-population split and radii injection behave as specified.")
            return 0
        return 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
