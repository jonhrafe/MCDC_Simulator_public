#!/usr/bin/env python3
"""
Permeable-membrane equilibrium (volume-preservation) regression test.

A single permeable ANALYTIC sphere (Sphere class, not a PLY mesh) sits in a periodic
cubic voxel whose side is chosen so the intracellular volume fraction is f = 0.5.
Intra and extra diffusivities differ by 10x (Di = 0.1, De = 1.0 um^2/ms), so the two
compartments use step lengths that differ by sqrt(10). Walkers are seeded uniformly.

For the Powles transit probabilities a uniform density is stationary, hence:
  1. the fraction of walkers inside the sphere must stay at f for the whole run
     (no drift),                                   -> |f(t) - f| <= TOL_F
  2. no walker may be discarded as an "illegal crossing" (legitimate crossings are
     not leaks)                                    -> eliminated <= MAX_ELIM_FRAC * N
Before the permeable-crossing fix (1-nm intra band in isInIntra + hit acceptance
1 nm beyond the step end) this setup drifts by several per cent within 20 ms and
discards thousands of walkers. See PERMEABLE_CROSSING_FIX.md.

`--substrate mesh` repeats the same check with a triangulated (PLY) sphere of the same
volume fraction, tracking the same per-step intra density through the mesh collision path.

Exit 0 = pass, 1 = fail. CPU only; ~2 min with 8 threads.
"""

import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

R_UM = 5.0
F_TARGET = 0.5
A_UM = (4.0 / 3.0 * math.pi * R_UM**3 / F_TARGET) ** (1.0 / 3.0)   # 10.1549 um
N = 40000
T = 40000                 # dt = 0.5 us
DURATION = 0.020          # s
TOL_F = 0.010             # 4 binomial SDs for N = 40000 (sd = 0.0025)
TOL_RELAX = 0.03          # max |f_mesh(t) - f_analytic(t)| for the all-intra relaxation (observed ~0.01)
MIN_RELAX = 0.15          # the all-intra start must have relaxed by at least this much by the end
MAX_ELIM_FRAC = 1e-3      # fixed build: 0 eliminated; unfixed build: ~2-4 %


def icosphere(radius, subdiv=3):
    """Outward-oriented icosphere (vertices, faces); returns also its enclosed volume."""
    t = (1.0 + math.sqrt(5.0)) / 2.0
    v = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
         (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    v = [tuple(c / math.sqrt(sum(x * x for x in p)) for c in p) for p in v]
    f = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4),
         (11, 10, 2), (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8),
         (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(subdiv):
        cache, nf = {}, []
        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = tuple((v[a][i] + v[b][i]) / 2 for i in range(3))
                n = math.sqrt(sum(x * x for x in m))
                v.append(tuple(x / n for x in m)); cache[k] = len(v) - 1
            return cache[k]
        for a, b, c in f:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nf += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        f = nf
    v = [tuple(radius * c for c in p) for p in v]
    vol = 0.0
    for i, (a, b, c) in enumerate(f):                      # enforce outward winding
        A, B, C = v[a], v[b], v[c]
        n = ((B[1]-A[1])*(C[2]-A[2]) - (B[2]-A[2])*(C[1]-A[1]),
             (B[2]-A[2])*(C[0]-A[0]) - (B[0]-A[0])*(C[2]-A[2]),
             (B[0]-A[0])*(C[1]-A[1]) - (B[1]-A[1])*(C[0]-A[0]))
        if sum(n[k] * (A[k] + B[k] + C[k]) for k in range(3)) < 0:
            f[i] = (a, c, b)
        a, b, c = f[i]; A, B, C = v[a], v[b], v[c]
        vol += (A[0]*(B[1]*C[2]-B[2]*C[1]) - A[1]*(B[0]*C[2]-B[2]*C[0]) + A[2]*(B[0]*C[1]-B[1]*C[0])) / 6.0
    return v, f, vol


def write_inputs(work, scheme_src, substrate="sphere", seeding="uniform"):
    side_um = A_UM
    if substrate == "mesh":
        verts, faces, vol = icosphere(R_UM)
        ply = os.path.join(work, "sphere.ply")
        with open(ply, "w") as fh:
            fh.write("ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\n"
                     "property float z\nelement face %d\nproperty list uchar int vertex_indices\n"
                     "end_header\n" % (len(verts), len(faces)))
            for p in verts:
                fh.write("%.9f %.9f %.9f\n" % p)
            for a, b, c in faces:
                fh.write("3 %d %d %d\n" % (a, b, c))
        side_um = (vol / F_TARGET) ** (1.0 / 3.0)           # f = 0.5 for the FACETED volume
        sph = os.path.join(work, "mesh.list")
        with open(sph, "w") as fh:
            # <ply> <scale m/unit> <d_intra m^2/s> <t2 s> <permeability m/s>
            fh.write("%s 1e-6 1.0e-10 1e9 5.0e-5\n" % ply)
        obstacle = "ply_extended_file_list %s" % sph
    else:
        sph = os.path.join(work, "sphere.txt")
        with open(sph, "w") as fh:
            fh.write("1e-6\n0 0 0 %.6f 0 1e9\n" % R_UM)  # scale(um->m); x y z r kappa(0=global) T2
        obstacle = "spheres_list %s" % sph
    scheme = os.path.join(work, "pgse.scheme")
    shutil.copy(scheme_src, scheme)
    h = 0.5 * side_um * 1e-6
    conf = os.path.join(work, "equil.conf")
    with open(conf, "w") as fh:
        fh.write(f"""N {N}
T {T}
duration {DURATION}
diff_intra 1.0e-10
diff_extra 1.0e-9
permeability 5.0e-5
exp_prefix {os.path.join(work, 'equil')}
scheme_file {scheme}
write_txt 1
write_bin 0
write_traj_file 0
deportation 1
discard_stucks 1
boundary_conditions periodic
<obstacle>
{obstacle}
</obstacle>
<voxels>
{-h:.9e} {-h:.9e} {-h:.9e}
{h:.9e} {h:.9e} {h:.9e}
</voxels>
{'ini_walkers_pos intra' if seeding == 'intra' else ''}
num_process 8
seed 20261021
<END>
""")
    return conf


def intra_fraction_series(work):
    files = sorted(glob.glob(os.path.join(work, "equil_*_particle_positions.txt")))
    if not files:
        raise RuntimeError("no *_particle_positions.txt written")
    total = None
    for f in files:
        vals = []
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                vals.append(int(line.split()[1]))
        total = vals if total is None else [a + b for a, b in zip(total, vals)]
    return [v / N for v in total[1:]]          # row 0 is the pre-seeding placeholder


def eliminated(work):
    txt = open(os.path.join(work, "equil_simulation_info.txt")).read()
    txt = re.sub(r"\x1b\[[0-9;]*m", "", txt)
    m = re.search(r"eliminated due crossings:\s*(\d+)", txt)
    return int(m.group(1)) if m else 0


def run_case(a, substrate, seeding):
    """Run one simulation in its own temp dir; return (intra-fraction series, eliminated)."""
    work = tempfile.mkdtemp(prefix="mcdc_equil_")
    try:
        conf = write_inputs(work, a.scheme, substrate, seeding)
        res = subprocess.run([a.bin, "--conf", conf], cwd=work,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if res.returncode != 0 or "[ERROR]" in (res.stdout or ""):
            sys.stderr.write(res.stdout or "")
            raise RuntimeError("simulator failed")
        return intra_fraction_series(work), eliminated(work)
    finally:
        if a.keep:
            print("work dir:", work)
        else:
            shutil.rmtree(work, ignore_errors=True)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bin", default=os.path.join(repo_root, "MC-DC_Simulator"))
    p.add_argument("--scheme", default=os.path.join(here, "accuracy", "PGSE_sample_scheme.scheme"))
    p.add_argument("--substrate", choices=("sphere", "mesh"), default="sphere",
                   help="analytic Sphere (default) or an equivalent triangulated PLY sphere")
    p.add_argument("--keep", action="store_true", help="keep the temporary work dirs")
    a = p.parse_args()
    ok = True
    tag = "equilibrium/" + a.substrate

    # --- 1. uniform seeding: the density must stay at f (no drift, no false discards)
    f, elim = run_case(a, a.substrate, "uniform")
    n10 = max(1, len(f) // 10)
    f_max_dev = max(abs(x - F_TARGET) for x in f)
    ok_f, ok_e = f_max_dev <= TOL_F, elim <= MAX_ELIM_FRAC * N
    print(f"[{tag}] uniform: f target={F_TARGET:.4f} first10%={sum(f[:n10])/n10:.4f} "
          f"last10%={sum(f[-n10:])/n10:.4f} max|f-f0|={f_max_dev:.4f} (tol {TOL_F})  {'OK' if ok_f else 'FAIL'}")
    print(f"[{tag}] uniform: eliminated as illegal crossings: {elim} of {N} "
          f"(max {int(MAX_ELIM_FRAC * N)})  {'OK' if ok_e else 'FAIL'}")
    ok &= ok_f and ok_e

    # --- 2. mesh only: an all-intra start must relax towards f exactly like the (validated)
    # analytic sphere. A uniform start is stationary even for an impermeable membrane, so this
    # is what proves the mesh actually exchanges, and at the right rate.
    if a.substrate == "mesh":
        fm, elim_m = run_case(a, "mesh", "intra")
        fa, elim_a = run_case(a, "sphere", "intra")
        dev = max(abs(x - y) for x, y in zip(fm, fa))
        relaxed = 1.0 - fm[-1]
        ok_r = dev <= TOL_RELAX and relaxed >= MIN_RELAX
        ok_e = elim_m <= MAX_ELIM_FRAC * N
        print(f"[{tag}] intra-seeded: f_mesh(end)={fm[-1]:.4f} f_analytic(end)={fa[-1]:.4f} "
              f"max|f_mesh-f_analytic|={dev:.4f} (tol {TOL_RELAX}), relaxed by {relaxed:.3f} "
              f"(min {MIN_RELAX})  {'OK' if ok_r else 'FAIL'}")
        print(f"[{tag}] intra-seeded: mesh eliminated: {elim_m} of {N}  {'OK' if ok_e else 'FAIL'}")
        ok &= ok_r and ok_e

    if ok:
        print("PASS: permeable compartments stay in equilibrium (volume preserved, no false discards).")
        return 0
    print("FAIL: permeable compartments not in equilibrium.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
