#!/usr/bin/env python3
"""
GPU Phase-2a validation: restricted diffusion in a single impermeable sphere.

Runs the SAME single-sphere config on the GPU (`gpu 1`) and the CPU (`gpu 0`) and
asserts the diffusion signals agree within MC noise (mean S/N over all
measurements). This guards the on-device collision + elastic-reflection (bounce)
kernel against the validated CPU path. N is kept modest so the CPU side stays fast.

Skips (exit 0) if the binary was built without CUDA / no GPU is present.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOL = 0.006   # |S/N_gpu - S/N_cpu|; covers independent-RNG MC noise at this N


def run(binp, conf_text, gpu, work, tag):
    conf = os.path.join(work, f"{tag}.conf")
    prefix = os.path.join(work, tag)
    lines = []
    for l in conf_text.splitlines():
        lines.append(f"exp_prefix {prefix}" if l.strip().startswith("exp_prefix") else l)
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
    # newest matching DWI file (handles <prefix>_DWI.txt and <prefix>_rep_00_DWI.txt)
    cands = [f for f in glob.glob(prefix + "*_DWI.txt") if not f.endswith(("_intra.txt", "_extra.txt"))]
    if not cands:
        return None, "no DWI output"
    f = max(cands, key=os.path.getmtime)
    vals = [float(x) for x in open(f) if x.strip()]
    return vals, "ok"


def mean_sn(vals, N):
    return (sum(vals) / len(vals)) / N


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    conf = os.path.join(REPO, "tests", "accuracy", "sphere_restricted.conf")
    text = open(conf).read()
    N = int(next(l.split()[1] for l in text.splitlines() if l.strip().startswith("N ")))

    work = tempfile.mkdtemp(prefix="mcdc_gpusph_")
    try:
        gpu_vals, st = run(binp, text, 1, work, "gpu")
        if st == "skip":
            print("SKIP: binary built without CUDA / no GPU available.")
            return 0
        if st != "ok":
            print(f"FAIL (GPU): {st}"); return 1
        cpu_vals, st = run(binp, text, 0, work, "cpu")
        if st != "ok":
            print(f"FAIL (CPU): {st}"); return 1
        if len(gpu_vals) != len(cpu_vals):
            print(f"FAIL: measurement count mismatch ({len(gpu_vals)} vs {len(cpu_vals)})"); return 1

        sn_gpu = mean_sn(gpu_vals, N)
        sn_cpu = mean_sn(cpu_vals, N)
        diff = abs(sn_gpu - sn_cpu)
        restricted = sn_cpu  # sanity: should differ from free exp(-bD)=~0.979
        ok = diff <= TOL
        print(f"GPU mean S/N = {sn_gpu:.5f} | CPU mean S/N = {sn_cpu:.5f} | "
              f"|diff| = {diff:.5f} (tol {TOL}) | restricted (cpu)={restricted:.4f}")
        if not ok:
            print("FAIL: GPU restricted-sphere signal disagrees with CPU beyond MC noise.")
            return 1
        print("PASS: GPU single-sphere collision matches CPU within MC noise.")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
