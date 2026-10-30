#!/usr/bin/env python3
"""
Substrate and run-control options, each checked against a known answer.

  U  SI vs use_mm_ms: the same substrate written in SI (m, s, T/m, list T2 in s, list scale in m per
     file unit) and in mm/ms (use_mm_ms 1, list scale in mm per file unit, T2 in ms, scheme in T/mm and
     ms) gives the same signals (total, intra, extra), for a spheres_list, a cylinders_list and a
     ply_file_list_scale_permeability (permeable spheres and mesh, per-sphere T2), and for start
     positions read from an ini_walkers_file.
  L  PLY lists: 'ply + ply_scale', 'ply_file_list' and 'ply_file_list_scale_permeability' with
     permeability 0 describe the same mesh and give the same signal; 'ply_file_list_scale_permeability'
     with kappa gives the same signal as 'ply_file_list' + obstacle_permeability kappa.
  B  boundary_conditions periodic / symmetric: free diffusion in a 5 um box (walkers cross the walls
     many times) gives exp(-b D) within Monte Carlo noise for both.
  O  oriented_cylinders_list: a 7-column (oriented) list is rejected with an explanation; a 6-column
     list read through oriented_cylinders_list gives the same signal as cylinders_list.
  X  max_sim_time: a run that cannot finish stops after the given wall-clock seconds; the signal is the
     sum over the walkers simulated so far, which simulation_info reports.
  V  verbatim (flag, no value): writes the _perm_counters.txt diagnostic for permeable obstacles;
     without it the file is not written.
  C  coll_sphere_dist: accepted with a warning that it has no effect.

  --gpu (GPU build) adds:
  G  every start mode on the GPU, by the number of walkers starting inside the obstacles (impermeable,
     no T2, so it is the intra signal at b = 0) for a sphere, a cylinder and a torus mesh:
       intra = N, extra = 0, delta inside = N, spawning box inside the obstacle = N, spawning box
       outside = 0, regular grid = the CPU count exactly, uniform = the volume fraction (4 SE),
     plus U and B on the GPU, and the rejection of mixed substrates.

Exit 0 = all pass.
"""
import argparse, os, re, shutil, subprocess, sys, tempfile, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_init_options import torus, write_ply  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCHEME = os.path.join(REPO, "docs", "scheme_files", "PGSE_sample_scheme.scheme")
GAMMA = 2.6751525e8
FAILS = []


def check(ok, msg):
    print("  [%s] %s" % ("OK " if ok else "FAIL", msg))
    if not ok: FAILS.append(msg)


def run(binary, work, name, body, expect_fail=False):
    prefix = os.path.join(work, name)
    for f in os.listdir(work):
        if f.startswith(name + "_"): os.remove(os.path.join(work, f))
    c = prefix + ".conf"
    open(c, "w").write(body.replace("@PREFIX@", prefix))
    t0 = time.time()
    r = subprocess.run([binary, "--conf", c], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    failed = r.returncode != 0 or "[ERROR]" in r.stdout
    if failed and not expect_fail:
        sys.stdout.write(r.stdout[-2000:]); raise RuntimeError(name + " failed")
    return prefix, failed, r.stdout, time.time() - t0


def dwi(prefix, part=""):
    return np.loadtxt(prefix + "_DWI" + part + ".txt")


def scheme_b():
    s = np.loadtxt(SCHEME, skiprows=1)
    return (GAMMA * s[:, 3] * s[:, 5]) ** 2 * (s[:, 4] - s[:, 5] / 3)   # s/m^2


def mm_scheme(path):
    s = np.loadtxt(SCHEME, skiprows=1)
    s[:, 3] *= 1e-3; s[:, 4:] *= 1e3                                    # T/m -> T/mm, s -> ms
    np.savetxt(path, s, header="VERSION: STEJSKALTANNER", comments="", fmt="%.10g")


def conf(prefix, N, T, D, dur, L, extra, scheme=SCHEME, nproc=2, mm=False, half=10, seed=99):
    """L = metres (SI) or mm (mm/ms) per micrometre."""
    return ("N %d\nT %d\nduration %g\ndiffusivity %g\nexp_prefix %s\nscheme_file %s\n%s"
            "write_txt 1\nwrite_bin 0\nseparate_signals\n<voxels>\n%g %g %g\n%g %g %g\n</voxels>\n%s\n"
            "num_process %d\nseed %d\n<END>\n"
            % (N, T, dur, D, prefix, scheme, "use_mm_ms 1\n" if mm else "",
               -half * L, -half * L, -half * L, half * L, half * L, half * L, extra, nproc, seed))


def same(a, b, tol=1e-4):
    """max |a-b| relative to max |a| (the b = 0 signal)."""
    return np.abs(a - b).max() / max(np.abs(a).max(), 1e-30) < tol, np.abs(a - b).max() / max(np.abs(a).max(), 1e-30)


def mesh_volume(P, F):
    A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    return abs((A * np.cross(B, C)).sum()) / 6


# ---------------------------------------------------------------- cases
def case_units(b, w, gpu=""):
    tag = "GPU " if gpu else ""
    print("U: %sSI vs use_mm_ms 1" % tag)
    P, F = torus(4.0, 1.5); P = P + np.array([-4.0, -4.0, 0.0]); write_ply(os.path.join(w, "torus_um.ply"), P, F)
    mm_scheme(os.path.join(w, "scheme_mm.scheme"))
    files = {"sph": ("5 5 5 3 1e-5 %s", "spheres_list"), "cyl": ("-5 5 0 2 0 %s", "cylinders_list")}
    for u, L, t2 in (("si", 1e-6, "0.05"), ("mm", 1e-3, "50")):
        for g, (row, _) in files.items():
            open(os.path.join(w, "%s_%s.txt" % (g, u)), "w").write("%g\n%s\n" % (L, row % t2))
        open(os.path.join(w, "ply_%s.txt" % u), "w").write("%s %g 2e-5\n" % (os.path.join(w, "torus_um.ply"), L))
        np.savetxt(os.path.join(w, "pos_%s.txt" % u), np.array([[1, 2, 3], [-4, 0, 2], [6, -6, -1]]) * L)
    geos = {"sph": "spheres_list {w}/sph_{u}.txt", "cyl": "cylinders_list {w}/cyl_{u}.txt",
            "ply": "ply_file_list_scale_permeability {w}/ply_{u}.txt", "file": None}
    for g, obst in geos.items():
        res = {}
        for u, L in (("si", 1e-6), ("mm", 1e-3)):
            mm = u == "mm"
            if obst: extra = "<obstacle>\n%s\n</obstacle>\n<spawning_area>\n%g %g %g\n%g %g %g\n</spawning_area>" % (
                obst.format(w=w, u=u), -8 * L, -8 * L, -8 * L, 8 * L, 8 * L, 8 * L)
            else: extra = "ini_walkers_file %s/pos_%s.txt" % (w, u)
            pre, _, _, _ = run(b, w, "u%s_%s_%s" % (gpu[:3], g, u), conf(
                "@PREFIX@", 4000, 500, 1.5e-6 if mm else 1.5e-9, 50 if mm else 0.05, L, extra + "\n" + gpu,
                scheme=os.path.join(w, "scheme_mm.scheme") if mm else SCHEME, mm=mm))
            res[u] = [dwi(pre, p) for p in ("", "_intra", "_extra")]
        oks = [same(a, c) for a, c in zip(res["si"], res["mm"])]
        check(all(o for o, _ in oks) and res["si"][0][0] > 0,
              "U %s%-4s: SI == mm/ms (max relative difference total %.1e, intra %.1e, extra %.1e)"
              % (tag, g, *[d for _, d in oks]))


def case_ply_lists(b, w):
    print("L: ply / ply_file_list / ply_file_list_scale_permeability")
    P, F = torus(5.0, 2.0); ply = os.path.join(w, "t5.ply"); write_ply(ply, P, F)
    open(os.path.join(w, "l1.txt"), "w").write("# mesh scale\n%s 1e-6\n" % ply)
    open(os.path.join(w, "l2.txt"), "w").write("%s 1e-6 0\n" % ply)
    open(os.path.join(w, "l3.txt"), "w").write("%s 1e-6 3e-5\n" % ply)
    variants = {"ply": "ply %s\nply_scale 1e-6" % ply, "list": "ply_file_list %s/l1.txt" % w,
                "list_perm0": "ply_file_list_scale_permeability %s/l2.txt" % w,
                "list_perm": "ply_file_list_scale_permeability %s/l3.txt" % w,
                "list_global": "ply_file_list %s/l1.txt" % w}
    S = {}
    for k, o in variants.items():
        pre, _, _, _ = run(b, w, "l_" + k, conf("@PREFIX@", 3000, 500, 1.5e-9, 0.05, 1e-6, "<obstacle>\n%s\n</obstacle>%s" % (
            o, "\nobstacle_permeability 3e-5" if k == "list_global" else "")))
        S[k] = dwi(pre, "_intra")
    for k in ("list", "list_perm0"):
        ok, d = same(S["ply"], S[k]); check(ok, "L %s == ply + ply_scale (max relative difference %.1e)" % (k, d))
    ok, d = same(S["list_perm"], S["list_global"])
    moved = not same(S["list_perm"], S["ply"], 1e-3)[0]
    check(ok and moved, "L per-mesh permeability 3e-5 == global obstacle_permeability 3e-5 (max rel. diff %.1e) "
          "and differs from the impermeable mesh" % d)


def free_z(S, D, N):
    E = np.exp(-scheme_b() * D)
    sig = np.sqrt(np.maximum((1 + E ** 4) / 2 - E ** 2, 1e-12) / N)
    z = (S / N - E) / np.maximum(sig, 1.0 / N)
    return np.abs(z).max(), z.mean()


def case_boundaries(b, w, gpu=""):
    tag = "GPU " if gpu else ""
    print("B: %sboundary_conditions periodic / symmetric (free diffusion, 5 um box)" % tag)
    N, D = 6000, 2e-9
    for bc in ("periodic", "symmetric"):
        pre, _, _, _ = run(b, w, "b%s_%s" % (gpu[:3], bc), conf("@PREFIX@", N, 1000, D, 0.05, 1e-6,
                           "boundary_conditions %s\n%s" % (bc, gpu), half=2.5))
        zmax, zmean = free_z(dwi(pre), D, N)
        check(zmax < 5 and abs(zmean) < 1.0, "B %s%-9s: signal vs exp(-bD) over %d measurements: max |z| %.1f, mean z %+.2f"
              % (tag, bc, len(scheme_b()), zmax, zmean))


def case_oriented(b, w):
    print("O: oriented_cylinders_list")
    open(os.path.join(w, "ocyl7.txt"), "w").write("1e-6\n0 0 -5 0 0 5 3\n")
    _, failed, out, _ = run(b, w, "o7", conf("@PREFIX@", 100, 10, 1.5e-9, 0.05, 1e-6,
                            "<obstacle>\noriented_cylinders_list %s/ocyl7.txt\n</obstacle>" % w), expect_fail=True)
    check(failed and "oriented" in out, "O 7-column list rejected with an explanation")
    open(os.path.join(w, "ocyl6.txt"), "w").write("1e-6\n0 0 0 3 0 1e9\n")
    S = {}
    for kw in ("cylinders_list", "oriented_cylinders_list"):
        pre, _, _, _ = run(b, w, "o_" + kw, conf("@PREFIX@", 2000, 200, 1.5e-9, 0.05, 1e-6,
                                                "<obstacle>\n%s %s/ocyl6.txt\n</obstacle>" % (kw, w)))
        S[kw] = dwi(pre)
    ok, d = same(S["cylinders_list"], S["oriented_cylinders_list"], 1e-9)
    check(ok, "O 6-column list: oriented_cylinders_list == cylinders_list")


def case_max_time(b, w):
    print("X: max_sim_time")
    pre, _, out, dt = run(b, w, "x", conf("@PREFIX@", 3000000, 1000, 2e-9, 0.05, 1e-6, "max_sim_time 2"))
    info = open(pre + "_simulation_info.txt").read()
    m = re.search(r"Number of simulated particles: (\d+)", info)
    n = int(m.group(1)) if m else -1
    s0 = dwi(pre)[0]
    check(dt < 15 and 0 < n < 3000000 and abs(s0 - n) < 0.5,
          "X stopped after %.1f s with %d of 3000000 walkers simulated; b = 0 signal %.0f = that number" % (dt, n, s0))


def case_verbatim(b, w):
    print("V: verbatim (flag) and coll_sphere_dist")
    open(os.path.join(w, "vsph.txt"), "w").write("1e-6\n0 0 0 4 2e-5 1e9\n")
    obst = "<obstacle>\nspheres_list %s/vsph.txt\n</obstacle>" % w
    pre, _, _, _ = run(b, w, "v_on", conf("@PREFIX@", 500, 100, 1.5e-9, 0.05, 1e-6, obst + "\nverbatim", nproc=1))
    pre2, _, _, _ = run(b, w, "v_off", conf("@PREFIX@", 500, 100, 1.5e-9, 0.05, 1e-6, obst, nproc=1))
    check(os.path.exists(pre + "_perm_counters.txt") and not os.path.exists(pre2 + "_perm_counters.txt"),
          "V _perm_counters.txt written with verbatim, not without")
    _, failed, out, _ = run(b, w, "c", conf("@PREFIX@", 100, 10, 1.5e-9, 0.05, 1e-6, "coll_sphere_dist 1e-7"))
    check(not failed and "coll_sphere_dist has no effect" in out, "C coll_sphere_dist accepted with a 'no effect' warning")


def intra0(pre):
    return dwi(pre, "_intra")[0]


def case_gpu_starts(cpu, gpu, w):
    print("G: start modes on the GPU (walkers starting inside the obstacle; impermeable, no T2)")
    P, F = torus(5.0, 2.0); write_ply(os.path.join(w, "gt.ply"), P, F)
    open(os.path.join(w, "gs.txt"), "w").write("1e-6\n0 0 0 4 0 1e9\n")
    open(os.path.join(w, "gc.txt"), "w").write("1e-6\n0 0 0 4 0 1e9\n")
    V = 20.0 ** 3
    geos = {"sphere": ("spheres_list %s/gs.txt" % w, 4 / 3 * np.pi * 64 / V, (0, 0, 0)),
            "cylinder": ("cylinders_list %s/gc.txt" % w, np.pi * 16 * 20 / V, (0, 0, 0)),
            "mesh": ("ply_file_list %s/gm.txt" % w, mesh_volume(P, F) / V, (5, 0, 0))}
    open(os.path.join(w, "gm.txt"), "w").write("%s/gt.ply 1e-6\n" % w)
    N = 20000
    for g, (obst, frac, c) in geos.items():
        box_in = "<spawning_area>\n%g %g %g\n%g %g %g\n</spawning_area>" % tuple(np.r_[np.array(c) - 0.8, np.array(c) + 0.8] * 1e-6)
        box_out = "<spawning_area>\n7e-6 7e-6 7e-6\n9e-6 9e-6 9e-6\n</spawning_area>"
        modes = {"intra": ("ini_walkers_pos intra", N), "extra": ("ini_walkers_pos extra", 0),
                 "delta": ("ini_walkers_pos delta\n<delta>\n%g %g %g\n</delta>" % tuple(np.array(c) * 1e-6 + 1e-7), N),
                 "box_inside": (box_in, N), "box_outside": (box_out, 0),
                 "regular": ("regular_sampling", None), "uniform": ("", frac)}
        for m, (ini, want) in modes.items():
            body = lambda dev: conf("@PREFIX@", N, 20, 1.5e-9, 0.05, 1e-6, "<obstacle>\n%s\n</obstacle>\n%s\n%s" % (obst, ini, dev))
            pg, _, _, _ = run(gpu, w, "g_%s_%s_gpu" % (g, m), body("gpu 1"))
            ng = intra0(pg)
            if m == "regular":
                pc, _, _, _ = run(cpu, w, "g_%s_%s_cpu" % (g, m), body(""))
                nc = intra0(pc)
                check(abs(ng - nc) < 0.5 and ng > 0, "G %-8s %-11s: GPU %d == CPU %d" % (g, m, ng, nc))
            elif m == "uniform":
                se = np.sqrt(frac * (1 - frac) / N)
                check(abs(ng / N - frac) < 4 * se, "G %-8s %-11s: GPU fraction %.4f vs volume fraction %.4f (4 SE %.4f)"
                      % (g, m, ng / N, frac, 4 * se))
            else:
                check(abs(ng - want) < 0.5, "G %-8s %-11s: GPU %d (expected %d)" % (g, m, ng, want))


def case_gpu_rejects(gpu, w):
    print("G: GPU rejects mixed substrates and warns on max_sim_time")
    open(os.path.join(w, "ms.txt"), "w").write("1e-6\n5 5 5 3 0 1e9\n")
    open(os.path.join(w, "mc.txt"), "w").write("1e-6\n-5 -5 0 2 0 1e9\n")
    _, failed, out, _ = run(gpu, w, "gx", conf("@PREFIX@", 100, 10, 1.5e-9, 0.05, 1e-6,
                            "<obstacle>\nspheres_list %s/ms.txt\ncylinders_list %s/mc.txt\n</obstacle>\ngpu 1" % (w, w)),
                            expect_fail=True)
    r = subprocess.run([gpu, "--conf", os.path.join(w, "gx.conf")], cwd=w, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    check(failed and r.returncode != 0, "G mixed spheres + cylinders with gpu 1: rejected, exit code %d" % r.returncode)
    _, _, out, _ = run(gpu, w, "gt_", conf("@PREFIX@", 100, 10, 1.5e-9, 0.05, 1e-6, "gpu 1\nmax_sim_time 5"))
    check("max_sim_time" in out and "GPU" in out, "G max_sim_time with gpu 1: warned that it is ignored")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=os.path.join(REPO, "MC-DC_Simulator"))
    ap.add_argument("--gpu", action="store_true", help="also run the GPU cases (CUDA build)")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    b = os.path.abspath(a.bin)
    work = tempfile.mkdtemp(prefix="mcdc_substrate_")
    try:
        case_units(b, work); case_ply_lists(b, work); case_boundaries(b, work); case_oriented(b, work)
        case_max_time(b, work); case_verbatim(b, work)
        if a.gpu:
            case_gpu_starts(b, b, work); case_units(b, work, "gpu 1"); case_boundaries(b, work, "gpu 1")
            case_gpu_rejects(b, work)
    finally:
        if a.keep: print("outputs kept in " + work)
        else: shutil.rmtree(work, ignore_errors=True)
    print("FAIL: %d check(s)" % len(FAILS) if FAILS else "PASS")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
