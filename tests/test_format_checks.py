#!/usr/bin/env python3
"""
Negative tests for the configuration-format fail-checks.

Feeds deliberately malformed .conf files to the simulator and asserts each one is
REJECTED (non-zero exit AND an "[ERROR]"/assertion message), while valid/benign
variants are ACCEPTED (exit 0). This guards the input validators (tag balance,
parameter ranges, PLY-list column counts, PLY property-vector alignment, and the
<bounding_box> obstacle-type restriction -- auto-voxel is only valid with
spheres/PLY, so cylinders or gamma packings must be rejected) and the P0.4 fix
that makes validation enforce in release builds (errors used to be emitted via
assert(), a no-op under -DNDEBUG, so bad configs slipped through).

Also: options whose value is missing or not a number (the parse used to stop silently), the
diffusivity checks (none given, explicit d_intra 0, t2_intra 0) and the d_intra / d_extra precedence
over `diffusivity` (read back from the run's parameter summary), an intra / extra start whose sampling
area does not reach the compartment (must exit with an error, not loop), and <bounding_box> without
obstacles.

Runs with a tiny N/T so accepted cases finish instantly. Exit 0 = all pass.
"""

import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "MC-DC_Simulator")
MESH = "tests/accuracy/meshes/unitMesh.ply"
SCHEME = "tests/accuracy/PGSE_sample_scheme.scheme"
CYL = "tests/accuracy/single_cylinder.list"

# A minimal VALID config (relative paths resolve from the repo root / workdir).
BASE = """N 4
T 4
duration 0.05
diffusivity 0.6e-9
exp_prefix {prefix}
scheme_file {scheme}
<obstacle>
ply {mesh}
ply_scale 1
</obstacle>
<voxels>
-1.1e-6 -1.1e-6 -1.1e-6
 1.1e-6  1.1e-6  1.1e-6
</voxels>
ini_walkers_pos intra
num_process 1
<END>
"""

# <bounding_box> auto-derives the voxel from obstacle extents, but is only
# supported for spheres / PLY meshes -- it is REJECTED (hard error) when CYLINDERS
# or GAMMA packings are present. These two configs must each be rejected.
BBOX_CYL = """N 4
T 4
duration 0.05
diffusivity 0.6e-9
exp_prefix {prefix}
scheme_file {scheme}
<obstacle>
cylinders_list {cyl}
</obstacle>
<bounding_box>
num_process 1
<END>
"""

BBOX_GAMMA = """N 4
T 4
duration 0.05
diffusivity 0.6e-9
exp_prefix {prefix}
scheme_file {scheme}
<obstacle>
<sphere_gamma_packing>
alpha 1.5
beta 0.5e-6
icvf 0.40
num_spheres 20
</sphere_gamma_packing>
</obstacle>
<bounding_box>
num_process 1
<END>
"""


def run(conf_text, work):
    conf = os.path.join(work, "t.conf")
    with open(conf, "w") as fh:
        fh.write(conf_text)
    res = subprocess.run([BIN, "--conf", conf], cwd=REPO,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                         timeout=120)
    out = res.stdout or ""
    err = ("[ERROR]" in out) or ("Assertion" in out)
    return res.returncode, err, out


def main():
    if not os.path.exists(BIN):
        print(f"FAIL: binary not found: {BIN}", file=sys.stderr)
        return 1

    work = tempfile.mkdtemp(prefix="mcdc_fmt_")
    try:
        base = BASE.format(prefix=os.path.join(work, "o"), scheme=SCHEME, mesh=MESH)

        # list files for the PLY-list cases
        ext4 = os.path.join(work, "ext4.list")          # extended list, 4 cols (missing permeability)
        with open(ext4, "w") as fh:
            fh.write(f"{MESH} 1 0.6e-9 0.08\n")
        fl1 = os.path.join(work, "fl1.list")            # simple list, 1 col (missing scale)
        with open(fl1, "w") as fh:
            fh.write(f"{MESH}\n")

        def obstacle_block(directive):
            return ("<obstacle>\n" + directive + "\n</obstacle>")

        # one sphere (r = 1 um at the origin, kappa 0, T2 0 = use t2_intra) for the placement cases
        sph = os.path.join(work, "one_sphere.list")
        with open(sph, "w") as fh:
            fh.write("1e-6\n0 0 0 1 0 0\n")
        sphere_base = base.replace(obstacle_block("ply " + MESH + "\nply_scale 1"),
                                   obstacle_block("<spheres_list>\n" + sph + "\n</spheres_list>"))
        # a spawning box in a corner of the voxel that the sphere does not reach (|corner| = 1.56 um)
        corner = "<spawning_area>\n0.9e-6 0.9e-6 0.9e-6\n1.1e-6 1.1e-6 1.1e-6\n</spawning_area>\n"
        D = "diffusivity 0.6e-9\n"

        # (name, conf_text, should_reject[, substring the output must contain])
        cases = [
            ("valid (control)",            base, False),
            ("comment with tags (benign)", base.replace("N 4", "# <obstacle> <delta> <voxels> note\nN 4"), False),
            ("unknown parameter (benign)", base.replace("N 4", "N 4\nfooparam 7"), False),
            ("missing </obstacle>",        base.replace("</obstacle>\n", ""), True),
            ("missing </voxels>",          base.replace("</voxels>\n", ""), True),
            ("zero N",                     base.replace("N 4", "N 0"), True),
            ("zero T",                     base.replace("T 4", "T 0"), True),
            ("inline ply, no ply_scale",   base.replace("ply_scale 1\n", ""), True),
            ("ply_extended_file_list 4cols",
                base.replace(obstacle_block("ply " + MESH + "\nply_scale 1"),
                             obstacle_block("ply_extended_file_list " + ext4)), True),
            ("ply_file_list missing scale",
                base.replace(obstacle_block("ply " + MESH + "\nply_scale 1"),
                             obstacle_block("ply_file_list " + fl1)), True),
            ("bounding_box + cylinder",
                BBOX_CYL.format(prefix=os.path.join(work, "o"), scheme=SCHEME, cyl=CYL), True),
            ("bounding_box + gamma packing",
                BBOX_GAMMA.format(prefix=os.path.join(work, "o"), scheme=SCHEME), True),
            # values: a missing or non-numeric value used to end the parse silently (the rest of
            # the .conf, voxel included, was ignored and the placement loop never returned)
            ("option without a value",     base.replace("<obstacle>", "deportation\n<obstacle>"), True),
            ("non-numeric value",          base.replace("N 4", "N four"), True),
            ("no diffusivity at all",      base.replace(D, ""), True),
            ("explicit d_intra 0",         base.replace(D, D + "d_intra 0\n"), True),
            ("t2_intra 0",                 base.replace(D, D + "t2_intra 0\n"), True),
            # d_intra / d_extra win over `diffusivity` whatever the order
            ("d_intra before diffusivity (precedence)",
                base.replace(D, "d_intra 1e-9\n" + D), False, "1.000000e-9 m^2/s", "0.600000e-9 m^2/s"),
            ("d_intra after diffusivity (precedence)",
                base.replace(D, D + "d_intra 1e-9\n"), False, "1.000000e-9 m^2/s", "0.600000e-9 m^2/s"),
            # placement: no position found -> error exit, not an endless loop (release builds)
            ("sphere list, T2 0 (control)", sphere_base, False),
            ("intra start, spawning area off the sphere",
                sphere_base.replace("ini_walkers_pos intra\n", "ini_walkers_pos intra\n" + corner), True,
                "Max. number of tries to find an intra-celular compartment"),
            ("extra start, voxel inside the mesh",
                base.replace("ini_walkers_pos intra", "ini_walkers_pos extra"), True,
                "Max. number of tries to find an extra-celular compartment"),
            # <bounding_box> needs spheres or meshes
            ("bounding_box without obstacles",
                base.replace(obstacle_block("ply " + MESH + "\nply_scale 1"), "<bounding_box>")
                    .replace("ini_walkers_pos intra", "ini_walkers_pos extra"), True,
                "No obstacles present"),
        ]

        ok = True
        print("Config-format fail-checks:")
        for name, text, should_reject, *expect in cases:
            rc, err, out = run(text, work)
            rejected = (rc != 0) and err
            passed = (rejected == should_reject) and all(e in out for e in expect)
            ok = ok and passed
            verdict = "OK  " if passed else "FAIL"
            kind = "reject" if should_reject else "accept"
            print(f"  {verdict} [{kind}] {name:32s} exit={rc} error={'y' if err else 'n'}")
            if not passed:
                # show a hint of why
                tail = " | ".join(l for l in out.splitlines() if "ERROR" in l or "Warning" in l)[:200]
                print(f"        -> {tail}")

        if ok:
            print("PASS: all format fail-checks behave as expected.")
            return 0
        print("FAIL: some format checks did not behave as expected.")
        return 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
