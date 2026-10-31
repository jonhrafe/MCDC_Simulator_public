#!/usr/bin/env python3
"""
Permeable equilibrium with steps comparable to the obstacle size.

A uniformly seeded permeable obstacle is in equilibrium: the fraction of walkers inside must stay
at the volume fraction f at every time. Here the obstacle radius is 1 um and the steps are
0.77 um (extra, D = 2e-9 m^2/s) and 0.39 um (intra, D = 0.5e-9 m^2/s), so a walker that crosses
the membrane can meet it again later in the same step. The engine must then cross or reflect it
again; finishing the step straight through after a crossing lets entering walkers run out the far
side (leaving ones cannot re-enter a convex obstacle) and drains the inside population (by 0.02 to
0.04 here). permeability_equilibrium cannot see this: its steps are ~1/50 of the radius.

For a sphere, a z-cylinder and a triangulated sphere (PLY), each in a periodic voxel, checks:
  - the inside fraction at 10 ms and at 50 ms equals f within 4 binomial SDs (with --gpu, also
    on the GPU engine when one is present);
  - no more than 0.1 % of the walkers are discarded.

Usage: test_permeability_coarse_step.py [--bin MC-DC_Simulator] [--gpu] [--substrate ...]
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from test_permeability_equilibrium import icosphere          # noqa: E402

R_UM = 1.0
SIDE_UM = 2.2                      # voxel side: f = 0.39 (sphere), 0.65 (cylinder)
D_INTRA, D_EXTRA = 0.5e-9, 2.0e-9  # m^2/s
KAPPA = 1.0e-4                     # m/s
DT = 5.0e-5                        # s: steps 0.77 um (extra), 0.39 um (intra)
DURATIONS = [0.010, 0.050]         # s
N = 100000
MAX_DROP = 1e-3                    # fraction of walkers allowed to be discarded


def write_ply(path, verts, faces):
    with open(path, "w") as fh:
        fh.write("ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\n"
                 "property float z\nelement face %d\nproperty list uchar int vertex_indices\n"
                 "end_header\n" % (len(verts), len(faces)))
        for p in verts:
            fh.write("%.9f %.9f %.9f\n" % p)
        for a, b, c in faces:
            fh.write("3 %d %d %d\n" % (a, b, c))


def write_inputs(work, substrate, duration, gpu):
    if substrate == "sphere":
        lst = os.path.join(work, "obs.list")
        with open(lst, "w") as fh:
            fh.write("1e-6\n0 0 0 %g %g 1e9\n" % (R_UM, KAPPA))       # x y z r kappa T2[s]
        obstacle = "spheres_list %s" % lst
        f = 4.0 / 3.0 * math.pi * R_UM**3 / SIDE_UM**3
    elif substrate == "cylinder":
        lst = os.path.join(work, "obs.list")
        with open(lst, "w") as fh:
            fh.write("1e-6\n0 0 0 %g %g 1e9\n" % (R_UM, KAPPA))       # along z
        obstacle = "cylinders_list %s" % lst
        f = math.pi * R_UM**2 / SIDE_UM**2
    else:
        verts, faces, vol = icosphere(R_UM)
        ply = os.path.join(work, "sphere.ply")
        write_ply(ply, verts, faces)
        lst = os.path.join(work, "mesh.list")
        with open(lst, "w") as fh:
            fh.write("%s 1e-6 %g 1e9 %g\n" % (ply, D_INTRA, KAPPA))   # ply scale d_intra T2 kappa
        obstacle = "ply_extended_file_list %s" % lst
        f = vol / SIDE_UM**3                                           # faceted volume
    scheme = os.path.join(work, "b0.scheme")
    with open(scheme, "w") as fh:
        fh.write("VERSION: STEJSKALTANNER\n1.000000 0.000000 0.000000 0.000000 0.001 0.0005 0.002\n")
    h = 0.5 * SIDE_UM * 1e-6
    conf = os.path.join(work, "eq.conf")
    with open(conf, "w") as fh:
        fh.write(f"""N {N}
T {int(round(duration / DT))}
duration {duration}
diff_intra {D_INTRA}
diff_extra {D_EXTRA}
t2_extra 1e9
exp_prefix {os.path.join(work, 'eq')}
scheme_file {scheme}
write_txt 1
write_bin 0
write_traj_file 0
separate_signals
boundary_conditions periodic
{"gpu 1" if gpu else ""}
<obstacle>
{obstacle}
</obstacle>
<voxels>
{-h:.9e} {-h:.9e} {-h:.9e}
{h:.9e} {h:.9e} {h:.9e}
</voxels>
num_process 8
seed 20261031
<END>
""")
    return conf, f


def run(binp, substrate, duration, gpu):
    """(inside fraction, walkers kept, f) | 'skip' | None."""
    work = tempfile.mkdtemp(prefix="mcdc_coarse_")
    try:
        conf, f = write_inputs(work, substrate, duration, gpu)
        res = subprocess.run([binp, "--conf", conf], cwd=work, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, timeout=1500)
        log = res.stdout or ""
        if gpu and ("WITHOUT CUDA" in log or "no CUDA-capable device" in log or "GPU path" not in log):
            return "skip"
        if res.returncode != 0 or "[ERROR]" in log:
            sys.stderr.write(log)
            return None
        try:
            si = float(open(os.path.join(work, "eq_DWI_intra.txt")).readline())
            se = float(open(os.path.join(work, "eq_DWI_extra.txt")).readline())
        except (OSError, ValueError) as e:
            sys.stderr.write("no intra/extra output: %s\n" % e)
            return None
        return si / (si + se), si + se, f
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bin", default=os.path.join(os.path.dirname(HERE), "MC-DC_Simulator"))
    p.add_argument("--gpu", action="store_true", help="also run the GPU engine (skipped without a GPU)")
    p.add_argument("--substrate", choices=("sphere", "cylinder", "mesh", "all"), default="all")
    a = p.parse_args()
    binp = os.path.abspath(a.bin)
    subs = ("sphere", "cylinder", "mesh") if a.substrate == "all" else (a.substrate,)
    engines = [("CPU", 0)] + ([("GPU", 1)] if a.gpu else [])

    fails = 0
    print("inside fraction vs volume fraction f (N = %d, tol = 4 binomial SDs)" % N)
    for sub in subs:
        for name, gpu in engines:
            for dur in DURATIONS:
                r = run(binp, sub, dur, gpu)
                if r == "skip":
                    print("  %-8s %s: SKIP (no GPU)" % (sub, name))
                    break
                if r is None:
                    print("  %-8s %s %4.0f ms: FAIL (simulation error)" % (sub, name, dur * 1e3))
                    fails += 1
                    continue
                frac, kept, f = r
                tol = 4.0 * math.sqrt(f * (1.0 - f) / N)
                drop = 1.0 - kept / N
                ok = abs(frac - f) <= tol and drop <= MAX_DROP
                fails += not ok
                print("  %-8s %s %4.0f ms: inside %.4f  f %.4f  |d| %.4f (tol %.4f)  dropped %d  %s"
                      % (sub, name, dur * 1e3, frac, f, abs(frac - f), tol, N - round(kept),
                         "PASS" if ok else "FAIL"))
    print("PASS" if fails == 0 else "FAIL: %d check(s)" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
