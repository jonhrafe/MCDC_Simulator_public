#!/usr/bin/env python3
"""
Logging / extra-output options: each output is checked against a known answer, in an SI .conf.

Case A, free diffusion (D = 0.6e-9 m^2/s, 50 ms, T = 200, 2 processes, periodic 100 um voxel), with
write_bin, write_traj_file, <log><positions>, <log><propagator>, log_phase_shift, log_opp, separate_signals:
  A1 write_bin      every .bfloat equals its .txt twin (float32)
  A2 traj header    .bhdr / .hdr.txt = (duration [ms], walkers in the file, number of logged times)
  A3 traj content   .traj == .traj.txt; positions logged only at the <positions> times; t=0 is the start;
                    per-axis displacement variance = 2 D t (internal units: mm, ms)
  A4 propagator     _propagator.txt row k = mean squared displacement along each direction at the k-th
                    logged time: equals the value recomputed from the trajectories, and 2 D t
  A5 phase shift    the |phase| histogram (3600 bins over [0, 2pi)) holds every walker and its mean
                    cos(phase) reproduces _DWI.txt / N
  A6 log_opp        the per-process _<p>_DWI.txt files add up to the merged _DWI.txt
  A7 separate       free diffusion: _DWI_intra = 0 and _DWI_extra = _DWI
Case B, full trajectory (no <positions>), 1 process: every step has length sqrt(6 D dt) exactly.
Case C, compute_volume: a 4 um sphere centred in a 10 um voxel, uniform start: the ICVF estimate in
  _simulation_info.txt matches (4/3 pi R^3)/V within binomial noise, and volume = ICVF * voxel volume.

Exit 0 = all pass.
"""
import argparse, glob, os, re, shutil, subprocess, sys, tempfile
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCHEME = os.path.join(REPO, "tests", "accuracy", "PGSE_sample_scheme.scheme")
D_SI, DUR_S, T = 0.6e-9, 0.050, 200   # duration covers the scheme's TE (45 ms)
D_INT, DUR_MS = D_SI * 1e3, DUR_S * 1e3          # mm^2/ms, ms
DIRS = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0]], float)
GAMMA = 2.6751525e8

FAILS = []


def check(ok, msg):
    print("  [%s] %s" % ("OK " if ok else "FAIL", msg))
    if not ok: FAILS.append(msg)


def run(binary, work, name, body):
    conf = os.path.join(work, name + ".conf")
    open(conf, "w").write(body)
    r = subprocess.run([binary, "--conf", conf], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0 or "[ERROR]" in r.stdout:
        sys.stdout.write(r.stdout[-2000:]); raise RuntimeError(name + " failed")
    return os.path.join(work, name)


def base(prefix, N, nproc, extra=""):
    return ("N %d\nT %d\nduration %g\ndiffusivity %g\nexp_prefix %s\nscheme_file %s\n%s"
            "num_process %d\nseed 12345\n<END>\n" % (N, T, DUR_S, D_SI, prefix, SCHEME, extra, nproc))


def bfloat(path): return np.fromfile(path, dtype="<f4")


def read_traj(prefix, p):
    hdr = bfloat("%s_%d.bhdr" % (prefix, p))
    n, nt = int(hdr[1]), int(hdr[2])
    pos = bfloat("%s_%d.traj" % (prefix, p))
    if pos.size == n * (nt + 1) * 3: nt += 1     # full trajectory: header holds T, the file T+1 positions
    pos = pos.reshape(n, nt, 3)
    return hdr, pos


def case_a(binary, work):
    print("Case A: free diffusion, all logging options, 2 processes")
    N, times = 4000, [0, 100, T]
    dirs = os.path.join(work, "dirs.txt"); np.savetxt(dirs, DIRS)
    body = base(os.path.join(work, "a"), N, 2,
                "write_txt 1\nwrite_bin 1\nwrite_traj_file 1\nseparate_signals\nlog_phase_shift\nlog_opp\n"
                "<log>\n<positions>\n0 100 t\n</positions>\n<propagator>\ndirections %s\n100 t\n</propagator>\n</log>\n"
                "<voxels>\n-50e-6 -50e-6 -50e-6\n50e-6 50e-6 50e-6\n</voxels>\n" % dirs)
    pre = run(binary, work, "a", body)

    # A1 bfloat == txt
    for suf in ("_DWI", "_DWI_intra", "_DWI_extra"):
        t, b = np.loadtxt(pre + suf + ".txt"), bfloat(pre + suf + ".bfloat")
        check(t.shape == b.shape and np.allclose(t, b, rtol=1e-5, atol=1e-6), "A1 %s.bfloat == %s.txt (txt has 6 significant digits)" % (suf, suf))
    ph_t = np.loadtxt(pre + "_phase_shift.txt"); ph_b = bfloat(pre + "_phase_shift.bfloat").reshape(ph_t.shape)
    check(np.array_equal(ph_t, ph_b), "A1 _phase_shift.bfloat == _phase_shift.txt (%dx%d)" % ph_t.shape)

    # A2/A3 trajectories
    pos_all, n_tot = [], 0
    for p in range(2):
        hdr, pos = read_traj(pre, p)
        th = np.loadtxt("%s_%d.hdr.txt" % (pre, p))
        check(np.allclose(hdr, th, rtol=1e-6), "A2 process %d: .bhdr == .hdr.txt = %s" % (p, th.tolist()))
        check(abs(hdr[0] - DUR_MS) < 1e-4 and int(hdr[2]) == len(times),
              "A2 process %d header: duration %.4g (expected %g ms), logged times %d (expected %d)"
              % (p, hdr[0], DUR_MS, int(hdr[2]), len(times)))
        txt = np.loadtxt("%s_%d.traj.txt" % (pre, p)).reshape(pos.shape)
        check(np.allclose(txt, pos, rtol=1e-5, atol=1e-9), "A3 process %d: .traj == .traj.txt" % p)
        pos_all.append(pos); n_tot += pos.shape[0]
    check(n_tot == N, "A2 walkers in the trajectory files: %d (expected %d)" % (n_tot, N))
    pos = np.concatenate(pos_all)
    disp = {k: pos[:, i, :] - pos[:, 0, :] for i, k in enumerate(times)}
    for k in (100, T):
        t_ms = DUR_MS * k / T
        var = disp[k].var(0) + disp[k].mean(0) ** 2
        rel = var / (2 * D_INT * t_ms) - 1
        tol = 4 * np.sqrt(2.0 / N)
        check(np.all(np.abs(rel) < tol), "A3 t=%d: <dx^2>/(2Dt) - 1 per axis = %s (tol %.3f)" % (k, np.round(rel, 4).tolist(), tol))
    ins = np.abs(pos[:, 0, :]) <= 0.05 + 1e-9
    check(ins.all(), "A3 every logged start position is inside the voxel (|x| <= 0.05 mm)")

    # A4 propagator
    prop = np.loadtxt(pre + "_propagator.txt", ndmin=2)
    check(prop.shape == (2, len(DIRS)), "A4 propagator shape %s (expected (2, %d))" % (prop.shape, len(DIRS)))
    pb = bfloat(pre + "_propagator.txt.bfloat")
    check(np.allclose(pb, prop.ravel(), rtol=1e-5), "A4 _propagator.txt.bfloat == _propagator.txt")
    u = DIRS / np.linalg.norm(DIRS, axis=1)[:, None]
    for row, k in enumerate((100, T)):
        msd = ((disp[k] @ u.T) ** 2).mean(0)
        check(np.allclose(prop[row], msd, rtol=2e-4), "A4 t=%d: propagator == MSD from the trajectories %s" % (k, np.round(msd / (2 * D_INT * DUR_MS * k / T), 4).tolist()))
        rel = prop[row] / (2 * D_INT * DUR_MS * k / T) - 1
        check(np.all(np.abs(rel) < 4 * np.sqrt(2.0 / N)), "A4 t=%d: propagator / (2Dt) - 1 = %s" % (k, np.round(rel, 4).tolist()))

    # A5 phase histogram vs signal
    dwi = np.loadtxt(pre + "_DWI.txt")
    check(np.all(ph_t.sum(1) == N), "A5 every histogram row holds all %d walkers" % N)
    centres = (np.arange(3600) + 0.5) * np.pi / 1800
    ecos = ph_t @ np.cos(centres) / N
    err = np.abs(ecos - dwi / N).max()
    check(err < 2e-3, "A5 mean cos(phase) from the histogram == DWI/N (max diff %.2e, bin width 0.1 deg)" % err)

    # A6 log_opp, A7 separate signals
    parts = sorted(glob.glob(pre + "_[0-9]_DWI.txt"))
    tot = sum(np.loadtxt(f) for f in parts)
    check(len(parts) == 2 and np.allclose(tot, dwi, rtol=1e-5, atol=1e-6), "A6 sum of %d per-process _DWI.txt == merged _DWI.txt" % len(parts))
    check(np.allclose(np.loadtxt(pre + "_DWI_intra.txt"), 0) and np.allclose(np.loadtxt(pre + "_DWI_extra.txt"), dwi),
          "A7 free diffusion: _DWI_intra = 0, _DWI_extra = _DWI")
    sch = np.loadtxt(SCHEME, skiprows=1)
    b = (GAMMA * sch[:, 3] * sch[:, 5]) ** 2 * (sch[:, 4] - sch[:, 5] / 3)
    E = np.exp(-b * D_SI)
    z = np.abs(dwi / N - E) / np.sqrt(np.clip(1 - E ** 2, 1e-6, None) / (2 * N))
    check(z.max() < 5, "A  free signal vs exp(-bD): max |z| = %.2f over %d rows (limit 5)" % (z.max(), len(z)))


def case_b(binary, work):
    print("Case B: full trajectory, step length")
    N = 200
    pre = run(binary, work, "b", base(os.path.join(work, "b"), N, 1,
              "write_txt 0\nwrite_bin 1\nwrite_traj_file 1\n<voxels>\n-50e-6 -50e-6 -50e-6\n50e-6 50e-6 50e-6\n</voxels>\n"))
    hdr, pos = read_traj(pre, 0)
    check(pos.shape == (N, T + 1, 3), "B .traj holds %s positions (expected (%d, %d, 3))" % (pos.shape, N, T + 1))
    step = np.linalg.norm(np.diff(pos, axis=1), axis=2)
    expect = np.sqrt(6 * D_INT * DUR_MS / T)
    check(np.allclose(step, expect, rtol=1e-4), "B every step = sqrt(6 D dt) = %.6g mm (min %.6g, max %.6g)" % (expect, step.min(), step.max()))


def case_c(binary, work):
    print("Case C: compute_volume, 4 um sphere in a 10 um voxel")
    N, R, L = 20000, 4e-6, 10e-6
    lst = os.path.join(work, "sphere.list"); open(lst, "w").write("1e-06\n0 0 0 4 0 1e9\n")   # x y z r permeability T2
    pre = run(binary, work, "c", base(os.path.join(work, "c"), N, 2,
              "write_txt 1\nwrite_bin 0\ncompute_volume\n<obstacle>\nspheres_list %s\n</obstacle>\n"
              "<voxels>\n-5e-6 -5e-6 -5e-6\n5e-6 5e-6 5e-6\n</voxels>\n" % lst))
    info = re.sub(r"\x1b\[[0-9;]*m", "", open(pre + "_simulation_info.txt").read())
    num = lambda key: float(re.search(key + r"[^\n]*\n\s*([-+0-9.eE]+)", info).group(1))
    icvf, vol, vtot = num("Estimated ICVF from sampling"), num(r"Estimated Intra-axonal volume"), num(r"Total sampling volume")
    f = 4 / 3 * np.pi * R ** 3 / L ** 3
    se = np.sqrt(f * (1 - f) / N)
    check(abs(icvf - f) < 4 * se + 1e-4, "C ICVF %.4f vs (4/3 pi R^3)/V = %.4f (4 SE = %.4f)" % (icvf, f, 4 * se))
    check(abs(vtot - (L * 1e3) ** 3) / (L * 1e3) ** 3 < 1e-5, "C total sampling volume %.4e mm^3 (expected %.4e)" % (vtot, (L * 1e3) ** 3))
    check(abs(vol - icvf * vtot) <= 1e-4 * vtot, "C intra volume %.4e mm^3 = ICVF x total" % vol)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=os.path.join(REPO, "MC-DC_Simulator"))
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    work = tempfile.mkdtemp(prefix="mcdc_log_")
    try:
        for c in (case_a, case_b, case_c):
            c(os.path.abspath(a.bin), work)
    finally:
        if a.keep: print("outputs kept in " + work)
        else: shutil.rmtree(work, ignore_errors=True)
    print("FAIL: %d check(s)" % len(FAILS) if FAILS else "PASS")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
