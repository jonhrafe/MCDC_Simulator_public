#!/usr/bin/env python3
"""
GPU free-diffusion validation.

Runs the simulator on the GPU (`gpu 1`) with no obstacles and asserts the PGSE
signal follows the analytic free-diffusion law S/S0 = exp(-b*D) for every
measurement, where b = (gamma*G*delta)^2 * (Delta - delta/3) is computed from the
scheme and D from the config (both SI). This guards the GPU dynamics+synthesis
path (the Phase-1 PoC) and stays valid as collision support is added later.

If the binary was built WITHOUT CUDA (or no GPU is present), the test SKIPS
(exit 0) rather than failing -- so CPU-only machines still pass the suite.

Exit 0 = pass or skip; non-zero = real failure.
"""
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GAMMA = 2.6751525e8   # rad/(s*T) -- matches giro=267.51525e3 rad/(ms*T)
RTOL = 0.025          # MC noise + finite-step bias. Fixed-length MC steps approximate
                      # Gaussian diffusion only to O(step), so S deviates from the continuum
                      # exp(-bD) by ~1/T, growing with b: ~1.7% at b=2000 (T=1000) on the
                      # modernized multi-shell scheme. Tightens to <0.3% at the old low b.


def read_conf_value(text, key):
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("#"):
            continue
        p = s.split()
        if p and p[0] == key:
            return p[1]
    return None


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    conf = os.path.join(REPO, "tests", "accuracy", "free_diffusion.conf")
    if not os.path.exists(binp):
        print(f"FAIL: binary not found: {binp}", file=sys.stderr); return 1
    text = open(conf).read()
    N = int(read_conf_value(text, "N"))
    D = float(read_conf_value(text, "diffusivity"))     # m^2/s
    scheme_rel = read_conf_value(text, "scheme_file")
    scheme_path = os.path.join(REPO, scheme_rel)

    work = tempfile.mkdtemp(prefix="mcdc_gpufree_")
    try:
        prefix = os.path.join(work, "free")
        conf_run = os.path.join(work, "run.conf")
        out = "\n".join((f"exp_prefix {prefix}" if l.strip().startswith("exp_prefix") else l)
                        for l in text.splitlines()) + "\n"
        open(conf_run, "w").write(out)

        res = subprocess.run([binp, "--conf", conf_run], cwd=REPO,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300)
        log = res.stdout or ""

        if "WITHOUT CUDA" in log or "CUDA error" in log or "no CUDA-capable device" in log:
            print("SKIP: binary built without CUDA / no GPU available.")
            return 0
        if res.returncode != 0 or "[ERROR]" in log:
            sys.stderr.write(log); print(f"FAIL: GPU run failed (exit {res.returncode})"); return 1

        dwi_files = glob.glob(prefix + "*DWI.txt")
        if not dwi_files:
            print("FAIL: no DWI output produced", file=sys.stderr); return 1
        dwi = [float(x) for x in open(dwi_files[0]) if x.strip()]

        # b per measurement from the scheme (all rows here share G/Delta/delta).
        rows = []
        for line in open(scheme_path):
            p = line.split()
            if len(p) >= 7:
                try:
                    rows.append([float(v) for v in p[:7]])
                except ValueError:
                    continue  # header
        if len(rows) != len(dwi):
            print(f"FAIL: scheme rows ({len(rows)}) != DWI entries ({len(dwi)})"); return 1

        worst = 0.0
        for i, (r, s) in enumerate(zip(rows, dwi)):
            G, Delta, delta = r[3], r[4], r[5]
            b = (GAMMA * G * delta) ** 2 * (Delta - delta / 3.0)
            expected = math.exp(-b * D)
            got = s / N
            rel = abs(got - expected) / expected
            worst = max(worst, rel)
            if rel > RTOL:
                print(f"FAIL: meas {i}: S/N={got:.5f} expected exp(-bD)={expected:.5f} "
                      f"(b={b:.4g}, rel={rel:.4f} > {RTOL})")
                return 1
        print(f"PASS: GPU free diffusion matches exp(-bD) for all {len(dwi)} measurements "
              f"(worst rel={worst:.4f}, exp(-bD)~{math.exp(-b*D):.4f}, N={N}).")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
