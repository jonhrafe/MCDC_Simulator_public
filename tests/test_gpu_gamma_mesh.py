#!/usr/bin/env python3
"""
GPU mesh path on a PACKED substrate: GPU-mesh vs CPU-mesh (intra-only).

Regression guard for the GPU triangle-mesh kernel on a dense, many-obstacle PLY.
It reproduces (in miniature) the bug where the mesh kernel

  * never honored `ini_walkers_pos intra` -- it seeded UNIFORMLY in the voxel, so
    an intra run started only the area-fraction of walkers inside the mesh and the
    rest free-diffused (the restricted signal collapsed), and
  * ran in float32 with an off-surface "nudge", leaking walkers through thin walls.

The test builds a gamma sphere packing (seed-deterministic), tessellates every
sphere -- including the periodic images the generator emits, so the mesh tiles the
periodic voxel -- into ONE PLY of icospheres, and runs it intra-only on the GPU
(`gpu 1`) and CPU (`gpu 0`). The CPU mesh path is the trusted reference; the GPU
mesh must match it within MC noise, and must retain (not leak away) the walkers.
Comparing GPU-vs-CPU on the *identical* mesh makes icosphere faceting cancel.

If the placement bug returns, the GPU intra signal collapses toward free diffusion
and diverges from the CPU by tens of pp. If wall-leak returns, walker retention
(b0) and the signal both drop on the GPU only.

Self-contained (writes its own scheme/configs, builds its own PLY to a temp dir).
Self-skips (exit 0) if the binary was built without CUDA / no GPU is present.
"""
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "MC-DC_Simulator")

# |mean S/N_gpu - mean S/N_cpu| tolerance on the identical mesh: independent-RNG MC
# noise at this N plus the float32(GPU)-vs-double(CPU) offset. The bug it guards was
# tens of pp, so 0.02 is comfortably tight. Matches test_gpu_restricted's TOL.
TOL = 0.02
# Max fraction of walkers the GPU may lose (flagged as leaked) relative to the CPU.
B0_TOL = 0.03

N, T, DUR, DIFF, SEED = 4000, 15000, 0.045, "2e-9", 777
GAMMA = ("<sphere_gamma_packing>\nalpha 3.0\nbeta 0.6e-6\nicvf 0.40\n"
         "num_spheres 60\nmin_radius 0.8e-6\n</sphere_gamma_packing>")


def write_scheme(path):
    # STEJSKALTANNER PGSE (Gx Gy Gz |G| Delta delta TE): b0 + x-ramp + z-ramp. Spheres
    # restrict isotropically, so both ramps probe the restriction (b up to ~2000 s/mm^2).
    lines = ["VERSION: STEJSKALTANNER", "0 0 0 0 0.025 0.010 0.045"]
    for G in (0.05, 0.09, 0.12):
        lines.append("1 0 0 %f 0.025 0.010 0.045" % G)
    for G in (0.05, 0.09, 0.12):
        lines.append("0 0 1 %f 0.025 0.010 0.045" % G)
    open(path, "w").write("\n".join(lines) + "\n")


def icosphere(subdiv):
    t = (1 + 5 ** 0.5) / 2
    V = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
         (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    V = [np.array(v) / np.linalg.norm(v) for v in V]
    F = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4),
         (11, 10, 2), (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8),
         (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    cache = {}

    def mid(a, b):
        k = (min(a, b), max(a, b))
        if k not in cache:
            m = V[a] + V[b]; m /= np.linalg.norm(m); V.append(m); cache[k] = len(V) - 1
        return cache[k]

    for _ in range(subdiv):
        F2 = []
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            F2 += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = F2
    return np.array(V), np.array(F, dtype=int)


def build_ply(list_path, ply_path, subdiv=2):
    sph = np.loadtxt(list_path, skiprows=1)
    sph = np.atleast_2d(sph)[:, :4]                 # x y z r (file units, = um here)
    uv, uf = icosphere(subdiv)
    VV, FF, off = [], [], 0
    for cx, cy, cz, r in sph:
        VV.append(uv * r + np.array([cx, cy, cz])); FF.append(uf + off); off += len(uv)
    V = np.vstack(VV); Fc = np.vstack(FF)
    with open(ply_path, "w") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write("element vertex %d\nproperty float x\nproperty float y\nproperty float z\n" % len(V))
        f.write("element face %d\nproperty list uchar int vertex_indices\nend_header\n" % len(Fc))
        for x, y, z in V:
            f.write("%.6f %.6f %.6f\n" % (x, y, z))
        for a, b, c in Fc:
            f.write("3 %d %d %d\n" % (a, b, c))
    return len(sph), len(Fc)


def parse_voxel_max_m(info_path):
    # simulation_info prints internal-unit (mm) voxel limits; convert mm -> m for the SI .conf.
    txt = open(info_path).read().splitlines()
    for i, l in enumerate(txt):
        if "Voxel limits" in l:
            vals = re.findall(r"\(([^)]+)\)", "\n".join(txt[i + 1:i + 3]))
            if len(vals) >= 2:
                mx = [float(x) for x in vals[1].split()]
                return mx[0] * 1e-3
    raise RuntimeError("voxel limits not found in " + info_path)


def run(conf_text, work, tag):
    conf = os.path.join(work, tag + ".conf")
    open(conf, "w").write(conf_text)
    res = subprocess.run([BIN, "--conf", conf], cwd=REPO, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, timeout=900)
    log = res.stdout or ""
    if "WITHOUT CUDA" in log or "CUDA error" in log or "no CUDA-capable device" in log:
        return None, "skip", log
    if res.returncode != 0 or "[ERROR]" in log:
        return None, "fail (exit %d)" % res.returncode, log
    prefix = os.path.join(work, tag)
    cands = [f for f in glob.glob(prefix + "*_DWI.txt") if not f.endswith(("_intra.txt", "_extra.txt"))]
    if not cands:
        return None, "no DWI", log
    vals = [float(x) for x in open(max(cands, key=os.path.getmtime)) if x.strip()]
    return vals, "ok", log


def mesh_conf(work, tag, gpu, ply, vox_m):
    prefix = os.path.join(work, tag)
    return ("N %d\nT %d\nduration %s\ndiffusivity %s\ngpu %d\nexp_prefix %s\n"
            "scheme_file %s\nwrite_txt 1\nwrite_bin 0\nt2_extra 1e9\n"
            "<obstacle>\nply %s\nply_scale 1e-6\n</obstacle>\nini_walkers_pos intra\n"
            "<voxels>\n0 0 0\n%g %g %g\n</voxels>\nnum_process 8\nseed %d\n<END>\n"
            % (N, T, DUR, DIFF, gpu, prefix, os.path.join(work, "perp.scheme"),
               ply, vox_m, vox_m, vox_m, SEED))


def main():
    work = tempfile.mkdtemp(prefix="mcdc_gmesh_")
    try:
        scheme = os.path.join(work, "perp.scheme"); write_scheme(scheme)

        # 1) Generate the packing (host-side, seed-deterministic) -> sphere list + voxel.
        gen_prefix = os.path.join(work, "gen")
        # N/T kept small but not degenerate (very tiny N*T aborts the substrate build);
        # this run only exists to emit the seed-deterministic packing list + voxel.
        gen = ("N 100\nT 100\nduration %s\ndiffusivity %s\ngpu 0\nexp_prefix %s\n"
               "scheme_file %s\nwrite_txt 1\nwrite_bin 0\nt2_extra 1e9\n"
               "<obstacle>\n%s\n</obstacle>\nini_walkers_pos intra\nnum_process 1\nseed %d\n<END>\n"
               % (DUR, DIFF, gen_prefix, scheme, GAMMA, SEED))
        _, st, log = run(gen, work, "gen")
        if st != "ok":
            sys.stderr.write(log); print("FAIL: packing generation: %s" % st); return 1
        lst = glob.glob(gen_prefix + "*_sphere_list.txt")
        if not lst:
            print("FAIL: no gamma sphere list produced"); return 1
        vox_m = parse_voxel_max_m(gen_prefix + "_simulation_info.txt")

        # 2) Tessellate the packing (incl. periodic images) into one PLY of icospheres.
        ply = os.path.join(work, "packing.ply")
        nsph, nfac = build_ply(lst[0], ply)
        print("packing: %d spheres (incl. periodic images) -> %d triangles; voxel %.2f um"
              % (nsph, nfac, vox_m * 1e6))

        # 3) GPU mesh first (so we can self-skip on non-CUDA builds before the CPU run).
        gpu_vals, st, log = run(mesh_conf(work, "gpu", 1, ply, vox_m), work, "gpu")
        if st == "skip":
            print("SKIP: binary built without CUDA / no GPU available."); return 0
        if st != "ok":
            sys.stderr.write(log); print("FAIL (GPU mesh): %s" % st); return 1

        # 4) CPU mesh reference (identical PLY, voxel, scheme, seed).
        cpu_vals, st, log = run(mesh_conf(work, "cpu", 0, ply, vox_m), work, "cpu")
        if st != "ok":
            sys.stderr.write(log); print("FAIL (CPU mesh): %s" % st); return 1

        if len(gpu_vals) != len(cpu_vals):
            print("FAIL: measurement count mismatch (%d vs %d)" % (len(gpu_vals), len(cpu_vals))); return 1

        sn_gpu = sum(gpu_vals) / len(gpu_vals) / N
        sn_cpu = sum(cpu_vals) / len(cpu_vals) / N
        diff = abs(sn_gpu - sn_cpu)
        b0_gpu, b0_cpu = gpu_vals[0] / N, cpu_vals[0] / N   # b0 = walker retention (S/N at b=0)
        last = cpu_vals[-1] / N                             # highest-b CPU signal
        print("GPU mesh S/N=%.5f  CPU mesh S/N=%.5f  |diff|=%.5f (tol %.3f)" % (sn_gpu, sn_cpu, diff, TOL))
        print("b0 retention: GPU=%.4f CPU=%.4f  |  CPU S/N at max b=%.4f" % (b0_gpu, b0_cpu, last))

        if any(math.isnan(v) for v in (sn_gpu, sn_cpu)):
            print("FAIL: NaN signal."); return 1
        if last > 0.97:
            print("FAIL: substrate is not restricting (no diffusion weighting) -- test is not exercising the mesh."); return 1
        if abs(b0_gpu - b0_cpu) > B0_TOL:
            print("FAIL: GPU mesh loses/leaks walkers vs CPU (b0 mismatch > %.0f%%)." % (B0_TOL * 100)); return 1
        if diff > TOL:
            print("FAIL: GPU mesh signal disagrees with CPU beyond MC noise (placement/leak regression?)."); return 1
        print("PASS: GPU mesh matches CPU mesh on a packed gamma substrate within MC noise.")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
