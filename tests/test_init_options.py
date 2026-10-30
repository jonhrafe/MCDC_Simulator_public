#!/usr/bin/env python3
"""
Walker-initialisation options, each checked against a known answer in an SI .conf. Start positions are
read from the trajectory file (write_traj_file 1, <log><positions> 0 t; positions in mm).

  D  <delta> + ini_walkers_pos delta: every walker starts at the point; at T the displacement
     variance is 2 D t (free diffusion), 2 processes.
  F  ini_walkers_file (positions in METRES in an SI .conf): walker g (global index over the processes)
     starts at file row g mod K, 2 processes.
  G  <regular_sampling>: N = 1000 on 2 processes fills the 10^3 cell centres of the voxel exactly
     once; N = 2000 (not a cube) spreads evenly over the 13^3 cells (balanced layers, centred mean).
  S  <spawning_area> with uniform / intra / extra starts among spheres: every start is inside the box,
     inside a sphere (intra) or outside all spheres (extra).
  M  intra / extra starts on a concave triangle mesh (a torus written by this test): every start agrees
     with a brute-force ray-parity inside test; a uniform start fills the mesh in proportion to its
     volume.
  P  the docs example walkerPositionsSpheres.conf (positions file among three spheres with different T2):
     the b = 0 intra / extra signals equal their exact values; also on the GPU in a GPU build.
  E  rejected configurations: ini_walkers_pos delta without <delta>; intra with <delta>; a file
     position outside the voxel.

Exit 0 = all pass.
"""
import argparse, os, shutil, struct, subprocess, sys, tempfile
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCHEME = os.path.join(REPO, "docs", "scheme_files", "PGSE_sample_scheme.scheme")
D_SI, DUR_S, T = 0.6e-9, 0.050, 100
FAILS = []


def check(ok, msg):
    print("  [%s] %s" % ("OK " if ok else "FAIL", msg))
    if not ok: FAILS.append(msg)


def conf(prefix, N, nproc, extra, voxel=10e-6, T_=T):
    return ("N %d\nT %d\nduration %g\ndiffusivity %g\nexp_prefix %s\nscheme_file %s\n"
            "write_txt 0\nwrite_bin 1\nwrite_traj_file 1\n<log>\n<positions>\n0 t\n</positions>\n</log>\n"
            "<voxels>\n%g %g %g\n%g %g %g\n</voxels>\n%s\nnum_process %d\nseed 12345\n<END>\n"
            % (N, T_, DUR_S, D_SI, prefix, SCHEME, -voxel, -voxel, -voxel, voxel, voxel, voxel, extra, nproc))


def run(binary, work, name, body, expect_fail=False):
    prefix = os.path.join(work, name)
    c = prefix + ".conf"
    open(c, "w").write(body.replace("@PREFIX@", prefix))
    r = subprocess.run([binary, "--conf", c], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    failed = r.returncode != 0 or "[ERROR]" in r.stdout
    if failed and not expect_fail:
        sys.stdout.write(r.stdout[-2000:]); raise RuntimeError(name + " failed")
    return prefix, failed, r.stdout


def starts(prefix, nproc):
    """(start, end) positions [mm] per process, from <prefix>_<p>.traj (2 logged times)."""
    out = []
    for p in range(nproc):
        hdr = np.fromfile("%s_%d.bhdr" % (prefix, p), "<f4")
        out.append(np.fromfile("%s_%d.traj" % (prefix, p), "<f4").reshape(int(hdr[1]), int(hdr[2]), 3))
    return out


# ---------------------------------------------------------------- torus mesh + parity ground truth
def torus(R, r, nu=48, nv=24):
    u = np.linspace(0, 2 * np.pi, nu, endpoint=False); v = np.linspace(0, 2 * np.pi, nv, endpoint=False)
    U, V = np.meshgrid(u, v, indexing="ij")
    P = np.stack([(R + r * np.cos(V)) * np.cos(U), (R + r * np.cos(V)) * np.sin(U), r * np.sin(V)], -1).reshape(-1, 3)
    idx = lambda i, j: (i % nu) * nv + (j % nv)
    F = []
    for i in range(nu):
        for j in range(nv):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            F += [(a, b, c), (a, c, d)]
    F = np.array(F)
    # orient every face outward (normal away from the tube centre line)
    A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    n = np.cross(B - A, C - A); cen = (A + B + C) / 3
    ring = np.stack([cen[:, 0], cen[:, 1], 0 * cen[:, 2]], -1)
    ring = R * ring / np.linalg.norm(ring, axis=1)[:, None]
    flip = ((cen - ring) * n).sum(1) < 0
    F[flip] = F[flip][:, [0, 2, 1]]
    return P, F


def write_ply(path, P, F):
    with open(path, "wb") as f:
        f.write(("ply\nformat binary_little_endian 1.0\nelement vertex %d\nproperty float x\nproperty float y\n"
                 "property float z\nelement face %d\nproperty list uchar int vertex_indices\nend_header\n" % (len(P), len(F))).encode())
        f.write(P.astype("<f4").tobytes())
        for t in F: f.write(struct.pack("<Biii", 3, *t))


def inside(Q, P, F):
    """brute-force ray parity (tilted ray) for points Q (k,3); same units as P."""
    d = np.array([0.0123, 0.0345, 1.0]); d /= np.linalg.norm(d)
    A, B, C = P[F[:, 0]].astype(float), P[F[:, 1]].astype(float), P[F[:, 2]].astype(float)
    e1, e2 = B - A, C - A
    pv = np.cross(d, e2); det = (e1 * pv).sum(1); inv = 1 / det
    res = np.zeros(len(Q), bool)
    for s in range(0, len(Q), 512):
        q = Q[s:s + 512, None, :]
        tv = q - A[None]
        u = (tv * pv[None]).sum(2) * inv
        qv = np.cross(tv, e1[None])
        v = (qv * d).sum(2) * inv
        t = (qv * e2[None]).sum(2) * inv
        hit = (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
        res[s:s + 512] = hit.sum(1) % 2 == 1
    return res


# ---------------------------------------------------------------- cases
def case_delta(b, w):
    print("D: <delta> start, 2 processes")
    N, pt = 600, np.array([1e-6, -2e-6, 0.5e-6])
    pre, _, _ = run(b, w, "d", conf("@PREFIX@", N, 2, "ini_walkers_pos delta\n<delta>\n%g %g %g\n</delta>" % tuple(pt)))
    S = np.concatenate(starts(pre, 2))
    check(S.shape[0] == N and np.allclose(S[:, 0], pt * 1e3, atol=1e-9), "D all %d walkers start at %s mm" % (len(S), (pt * 1e3).tolist()))
    msd = ((S[:, 1] - S[:, 0]) ** 2).mean(0) / (2 * D_SI * 1e3 * DUR_S * 1e3) - 1
    check(np.all(np.abs(msd) < 4 * np.sqrt(2 / N)), "D displacement variance / (2Dt) - 1 = %s" % np.round(msd, 3).tolist())


def case_file(b, w):
    print("F: ini_walkers_file (metres), 2 processes")
    K = np.array([[1e-6, 2e-6, 3e-6], [-4e-6, 5e-6, -6e-6], [7e-6, -8e-6, 9e-6]])
    f = os.path.join(w, "ini.txt"); np.savetxt(f, K)
    N = 7
    pre, _, _ = run(b, w, "f", conf("@PREFIX@", N, 2, "ini_walkers_file %s" % f))
    per = starts(pre, 2)
    S = np.concatenate([p[:, 0] for p in per])
    expect = K[np.arange(N) % len(K)] * 1e3
    check(np.allclose(S, expect, atol=1e-9), "F walker g starts at file row g mod %d, in mm (= metres x 1000), across processes" % len(K))


def case_regular(b, w):
    print("G: <regular_sampling>")
    pre, _, _ = run(b, w, "g1", conf("@PREFIX@", 1000, 2, "<regular_sampling>"))
    S = np.concatenate([p[:, 0] for p in starts(pre, 2)])
    centres = -0.01 + (np.arange(10) + 0.5) * 0.002
    grid = np.stack(np.meshgrid(centres, centres, centres, indexing="ij"), -1).reshape(-1, 3)
    near = np.argmin(((S[:, None, :] - grid[None]) ** 2).sum(2), axis=1)        # nearest cell centre
    err = np.abs(S - grid[near]).max()                                           # float32 positions
    check(len(np.unique(near)) == 1000 and err < 1e-8,
          "G N=1000, 2 processes: the 1000 starts are the 10^3 cell centres, each used once (max dev %.1e mm)" % err)
    pre, _, _ = run(b, w, "g2", conf("@PREFIX@", 2000, 1, "<regular_sampling>"))
    S = starts(pre, 1)[0][:, 0]
    layers = [np.unique(np.round(S[:, a], 7), return_counts=True)[1] for a in range(3)]
    bal = all(len(c) == 13 and c.max() - c.min() <= 0.1 * c.mean() for c in layers)
    check(bal, "G N=2000: 13 layers per axis, counts %s .. %s (balanced)" % ([c.min() for c in layers], [c.max() for c in layers]))
    check(np.all(np.abs(S.mean(0)) < 1e-3 * 0.02), "G N=2000: mean start position %s mm (centred)" % np.round(S.mean(0), 6).tolist())


def case_spawn(b, w):
    print("S: <spawning_area> with uniform / intra / extra starts among spheres")
    sph = np.array([[-4e-6, -4e-6, 0, 3e-6], [4e-6, 3e-6, 1e-6, 2.5e-6], [0, 4e-6, -5e-6, 2e-6]])
    lst = os.path.join(w, "spheres.list")
    with open(lst, "w") as f:
        f.write("1e-06\n")
        for s in sph: f.write("%g %g %g %g 0 1e9\n" % tuple(s * 1e6))
    box = np.array([[-6e-6, -6e-6, -6e-6], [5e-6, 6e-6, 4e-6]])
    spawn = "<spawning_area>\n%g %g %g\n%g %g %g\n</spawning_area>" % tuple(box.ravel())
    obst = "<obstacle>\nspheres_list %s\n</obstacle>" % lst
    for mode in ("uniform", "intra", "extra"):
        ini = "" if mode == "uniform" else "ini_walkers_pos %s" % mode
        pre, _, _ = run(b, w, "s_" + mode, conf("@PREFIX@", 1000, 2, "%s\n%s\n%s" % (obst, spawn, ini)))
        S = np.concatenate([p[:, 0] for p in starts(pre, 2)]) * 1e-3           # mm -> m
        inbox = np.all((S >= box[0] - 1e-12) & (S <= box[1] + 1e-12), 1)
        dist = np.linalg.norm(S[:, None, :] - sph[None, :, :3], axis=2) - sph[None, :, 3]
        intra = (dist < 0).any(1)
        msg = "S %-7s: %d/%d starts inside the spawning box" % (mode, inbox.sum(), len(S))
        if mode == "intra": ok = inbox.all() and intra.all(); msg += ", %d/%d inside a sphere" % (intra.sum(), len(S))
        elif mode == "extra": ok = inbox.all() and not intra.any(); msg += ", %d/%d outside all spheres" % ((~intra).sum(), len(S))
        else: ok = inbox.all(); msg += ", intra fraction %.3f" % intra.mean()
        check(ok, msg)


def case_mesh(b, w):
    print("M: intra / extra / uniform starts on a concave mesh (torus, R = 6 um, r = 2.5 um)")
    P, F = torus(6e-6, 2.5e-6)
    ply = os.path.join(w, "torus.ply"); write_ply(ply, P, F)
    obst = "<obstacle>\nply %s\nply_scale 1\n</obstacle>" % ply
    for mode in ("intra", "extra"):
        pre, _, _ = run(b, w, "m_" + mode, conf("@PREFIX@", 2000, 2, "%s\nini_walkers_pos %s" % (obst, mode)))
        S = np.concatenate([p[:, 0] for p in starts(pre, 2)]) * 1e-3
        ins = inside(S, P, F)
        good = ins if mode == "intra" else ~ins
        check(good.all(), "M %s: %d/%d starts on the right side of the mesh (ray parity)" % (mode, good.sum(), len(S)))
    N = 4000
    pre, _, _ = run(b, w, "m_uniform", conf("@PREFIX@", N, 2, obst))
    S = np.concatenate([p[:, 0] for p in starts(pre, 2)]) * 1e-3
    frac = inside(S, P, F).mean()
    rng = np.random.default_rng(1); Q = rng.uniform(-10e-6, 10e-6, (40000, 3))
    ref = inside(Q, P, F).mean()
    se = np.sqrt(ref * (1 - ref) / N)
    check(abs(frac - ref) < 4 * se, "M uniform: fraction of starts inside %.4f vs mesh volume fraction %.4f (4 SE %.4f)" % (frac, ref, 4 * se))



def case_file_example(b, w):
    """Runs docs/conf_file_examples/walkerPositionsSpheres.conf: 4 file rows (3 sphere centres with
    T2 0.05/0.10/0.20 s, 1 point outside), N = 1001 -> 251/250/250/250 walkers per row. Impermeable
    spheres, so the b = 0 intra and extra signals are known exactly. Also with gpu 1 in a GPU build."""
    print("P: positions file among spheres (docs example), CPU and GPU")
    import math
    want_i = 251 * math.exp(-0.045 / 0.05) + 250 * math.exp(-0.045 / 0.10) + 250 * math.exp(-0.045 / 0.20)
    ex = open(os.path.join(REPO, "docs", "conf_file_examples", "walkerPositionsSpheres.conf")).read()
    for dev in ("CPU", "GPU"):
        name = "p_" + dev.lower()
        body = ex.replace("exp_prefix instructions/demos/output/walker_positions_spheres", "exp_prefix @PREFIX@")
        if dev == "GPU": body = body.replace("# gpu 1", "gpu 1")
        prefix = os.path.join(w, name)
        c = prefix + ".conf"; open(c, "w").write(body.replace("@PREFIX@", prefix))
        r = subprocess.run([b, "--conf", c], cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if dev == "GPU" and "GPU path" not in r.stdout:
            print("  [SKIP] GPU: not a GPU build"); continue
        if r.returncode != 0 or "[ERROR]" in r.stdout:
            sys.stdout.write(r.stdout[-1500:]); check(False, "P %s: example run failed" % dev); continue
        si = float(np.loadtxt(prefix + "_DWI_intra.txt")[0]); se = float(np.loadtxt(prefix + "_DWI_extra.txt")[0])
        check(abs(si - want_i) < 1e-3 * want_i and abs(se - 250) < 1e-3,
              "P %s: b=0 intra %.3f (expected %.3f), extra %.3f (expected 250)" % (dev, si, want_i, se))

def case_errors(b, w):
    print("E: rejected configurations")
    _, failed, out = run(b, w, "e1", conf("@PREFIX@", 10, 1, "ini_walkers_pos delta"), expect_fail=True)
    check(failed and "<delta>" in out, "E ini_walkers_pos delta without <delta> is rejected")
    _, failed, _ = run(b, w, "e2", conf("@PREFIX@", 10, 1, "ini_walkers_pos intra\n<delta>\n0 0 0\n</delta>"), expect_fail=True)
    check(failed, "E ini_walkers_pos intra together with <delta> is rejected")
    _, failed, _ = run(b, w, "e4", conf("@PREFIX@", 10, 1, "<regular_sampling>\n<delta>\n0 0 0\n</delta>"), expect_fail=True)
    check(failed, "E <regular_sampling> together with <delta> is rejected")
    f0 = os.path.join(w, "inside.txt"); np.savetxt(f0, [[0, 0, 0]])
    _, failed, _ = run(b, w, "e5", conf("@PREFIX@", 10, 1, "ini_walkers_file %s\n<delta>\n0 0 0\n</delta>" % f0), expect_fail=True)
    check(failed, "E ini_walkers_file together with <delta> is rejected")
    f = os.path.join(w, "outside.txt"); np.savetxt(f, [[0, 0, 0], [20e-6, 0, 0]])
    _, failed, _ = run(b, w, "e3", conf("@PREFIX@", 10, 1, "ini_walkers_file %s" % f), expect_fail=True)
    check(failed, "E an ini_walkers_file position outside the voxel (20 um, voxel +-10 um) is rejected")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=os.path.join(REPO, "MC-DC_Simulator"))
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    work = tempfile.mkdtemp(prefix="mcdc_init_")
    try:
        for c in (case_delta, case_file, case_file_example, case_regular, case_spawn, case_mesh, case_errors):
            c(os.path.abspath(a.bin), work)
    finally:
        if a.keep: print("outputs kept in " + work)
        else: shutil.rmtree(work, ignore_errors=True)
    print("FAIL: %d check(s)" % len(FAILS) if FAILS else "PASS")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
