#!/usr/bin/env python3
"""
Negative tests for the configuration-format fail-checks.

Feeds deliberately malformed .conf files to the simulator and asserts each one is
REJECTED (non-zero exit AND an "[ERROR]"/assertion message), while valid/benign
variants are ACCEPTED (exit 0). This guards the input validators (tag balance,
parameter ranges, PLY-list column counts, PLY property-vector alignment) and the
P0.4 fix that makes validation enforce in release builds (errors used to be
emitted via assert(), a no-op under -DNDEBUG, so bad configs slipped through).

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

        # (name, conf_text, should_reject)
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
        ]

        ok = True
        print("Config-format fail-checks:")
        for name, text, should_reject in cases:
            rc, err, out = run(text, work)
            rejected = (rc != 0) and err
            passed = (rejected == should_reject)
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
