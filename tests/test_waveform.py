#!/usr/bin/env python3
"""
General-waveform protocol validation (GradientWaveform sequence).

Encodes a PGSE as an arbitrary gradient waveform (bipolar: +G for delta, then -G
for delta after Delta -- exactly the two-lobe Stejskal-Tanner shape) and runs it
on free diffusion, alongside the equivalent NATIVE PGSE. Both are normalized by
their own b0 (G=0) acquisition. Asserts the waveform signal matches BOTH the
analytic exp(-bD) and the native PGSE within MC + discretization noise. Guards the
WAVEFORM scheme path (scheme.type==WAVEFORM -> GradientWaveform) end to end. CPU-only.

Usage: test_waveform.py
"""
import glob
import math
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# MC noise at N=30k PLUS waveform discretization vs the analytic PGSE lobes (~0.3% at
# 1801 bins). Comfortably tight enough to catch a broken waveform encoding.
TOL = 0.02

# The encoded acquisition: b = (gamma*G*delta)^2 * (Delta - delta/3), free D.
# Must match tests/accuracy/{pgse_single,waveform_pgse}.scheme and *_free.conf.
G, DELTA, DELTA_SMALL, D = 0.05, 0.025, 0.005, 2.0e-9
GAMMA = 2.6751525e8  # rad/(s*T)
B = (GAMMA * G * DELTA_SMALL) ** 2 * (DELTA - DELTA_SMALL / 3.0)
EXPECTED = math.exp(-B * D)


def run(name):
    conf = os.path.join(REPO, "tests", "accuracy", name + ".conf")
    pref = None
    for l in open(conf):
        if l.strip().startswith("exp_prefix"):
            pref = os.path.join(REPO, l.split(None, 1)[1].strip())
    for old in glob.glob(pref + "*"):  # avoid reading a stale _rep_NN_ output
        try:
            os.remove(old)
        except OSError:
            pass
    res = subprocess.run([os.path.join(REPO, "MC-DC_Simulator"), "--conf", conf], cwd=REPO,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    if res.returncode != 0 or "[ERROR]" in (res.stdout or ""):
        sys.stderr.write(res.stdout or "")
        return None
    cands = [f for f in glob.glob(pref + "*_DWI.txt") if not f.endswith(("_intra.txt", "_extra.txt"))]
    if not cands:
        return None
    return [float(x) for x in open(max(cands, key=os.path.getmtime)) if x.strip()]


def main():
    norm = {}
    for name in ["pgse_free", "waveform_free"]:
        v = run(name)
        if v is None:
            print(f"FAIL: {name} run failed")
            return 1
        if len(v) < 2 or v[0] <= 0:
            print(f"FAIL: {name} expected b0+encoded (got {v})")
            return 1
        norm[name] = v[1] / v[0]
    pgse, wave = norm["pgse_free"], norm["waveform_free"]
    print(f"PGSE norm={pgse:.4f}  WAVEFORM norm={wave:.4f}  exp(-bD)={EXPECTED:.4f} "
          f"(b={B*1e-6:.1f} s/mm^2)")
    if pgse != pgse or wave != wave:
        print("FAIL: NaN signal")
        return 1
    if abs(wave - EXPECTED) > TOL:
        print(f"FAIL: waveform vs exp(-bD) |{wave-EXPECTED:.4f}| > {TOL}")
        return 1
    if abs(wave - pgse) > TOL:
        print(f"FAIL: waveform vs PGSE |{wave-pgse:.4f}| > {TOL}")
        return 1
    print("PASS: general-waveform protocol matches PGSE and exp(-bD).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
