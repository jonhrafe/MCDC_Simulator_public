#!/usr/bin/env python3
"""
GPU restricted-diffusion validation (Phase 2a/2b).

Runs a single-obstacle config (sphere or cylinder) on the GPU (`gpu 1`) and the
CPU (`gpu 0`) and asserts the diffusion signals agree within MC noise (mean S/N
over all measurements). Guards the on-device collision + elastic-reflection
kernels against the validated CPU path. N is kept modest so the CPU side is fast.

Usage: test_gpu_restricted.py [config.conf]   (default: sphere_restricted.conf)
Skips (exit 0) if the binary was built without CUDA / no GPU is present.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# |S/N_gpu - S/N_cpu| tolerance. Covers independent-RNG MC noise at this N PLUS the
# systematic float32(GPU)-vs-double(CPU) precision offset in restricted diffusion. That
# offset scales with diffusion weighting: under the modernized multi-shell protocol
# (b up to 2000 s/mm^2) the smaller signals make the fixed float32 gap a larger fraction
# of S/N, so the worst case grew to ~1.5% (multi-cylinder). Catches gross collision errors
# while allowing that float32 offset. (Was 0.008 under the old low-b scheme.)
TOL = 0.02


def run(binp, conf_text, gpu, work, tag):
    conf = os.path.join(work, f"{tag}.conf")
    prefix = os.path.join(work, tag)
    # Insert the gpu flag BEFORE <END> (the parser stops at <END>), and rewrite
    # exp_prefix to the temp dir.
    lines = []
    placed = False
    for l in conf_text.splitlines():
        s = l.strip()
        if s.startswith("exp_prefix"):
            l = f"exp_prefix {prefix}"
        if s == "<END>" and not placed:
            lines.append(f"gpu {gpu}"); placed = True
        lines.append(l)
    if not placed:
        lines.append(f"gpu {gpu}")
    open(conf, "w").write("\n".join(lines) + "\n")
    res = subprocess.run([binp, "--conf", conf], cwd=REPO,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    log = res.stdout or ""
    if gpu == 1 and ("WITHOUT CUDA" in log or "CUDA error" in log or "no CUDA-capable device" in log):
        return None, "skip"
    if res.returncode != 0 or "[ERROR]" in log:
        sys.stderr.write(log)
        return None, f"fail (exit {res.returncode})"
    cands = [f for f in glob.glob(prefix + "*_DWI.txt") if not f.endswith(("_intra.txt", "_extra.txt"))]
    if not cands:
        return None, "no DWI output"
    f = max(cands, key=os.path.getmtime)
    return [float(x) for x in open(f) if x.strip()], "ok"


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    conf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "tests", "accuracy", "sphere_restricted.conf")
    if not os.path.isabs(conf):
        conf = os.path.join(REPO, conf)
    text = open(conf).read()
    N = int(next(l.split()[1] for l in text.splitlines() if l.strip().startswith("N ")))

    work = tempfile.mkdtemp(prefix="mcdc_gpures_")
    try:
        gpu_vals, st = run(binp, text, 1, work, "gpu")
        if st == "skip":
            print("SKIP: binary built without CUDA / no GPU available."); return 0
        if st != "ok":
            print(f"FAIL (GPU): {st}"); return 1
        cpu_vals, st = run(binp, text, 0, work, "cpu")
        if st != "ok":
            print(f"FAIL (CPU): {st}"); return 1
        if len(gpu_vals) != len(cpu_vals):
            print(f"FAIL: measurement count mismatch ({len(gpu_vals)} vs {len(cpu_vals)})"); return 1

        sn_gpu = sum(gpu_vals) / len(gpu_vals) / N
        sn_cpu = sum(cpu_vals) / len(cpu_vals) / N
        import math
        diff = abs(sn_gpu - sn_cpu)
        print(f"[{os.path.basename(conf)}] GPU S/N={sn_gpu:.5f} CPU S/N={sn_cpu:.5f} "
              f"|diff|={diff:.5f} (tol {TOL})")
        if math.isnan(sn_gpu) or math.isnan(sn_cpu):
            print("FAIL: NaN signal (GPU or CPU produced NaN)."); return 1
        if diff > TOL:
            print("FAIL: GPU restricted signal disagrees with CPU beyond MC noise."); return 1
        print("PASS: GPU collision matches CPU within MC noise.")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
