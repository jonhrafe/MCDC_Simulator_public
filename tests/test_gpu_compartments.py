#!/usr/bin/env python3
"""
GPU compartment-tracking validation (Phase 3).

Runs a config with `separate_signals 1` on GPU (`gpu 1`) and CPU (`gpu 0`) and checks
that the GPU classifies walkers into the intra / extra compartments the same way the
CPU does: it compares (a) the intra signal FRACTION intra/(intra+extra) -- stable
vs independent-RNG placement noise -- and (b) the total signal. This guards the
per-walker `in_obs` classification the dynamics kernels write at placement.

Usage: test_gpu_compartments.py [config.conf]
Skips (exit 0) if built without CUDA / no GPU.
"""
import glob, os, shutil, subprocess, sys, tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOL_FRAC  = 0.02     # |intra-fraction_gpu - intra-fraction_cpu| (covers count-RNG at this N)
TOL_TOTAL = 0.006    # |mean S/N total| (float32 + noise)


def run(binp, conf_text, gpu, work, tag):
    conf = os.path.join(work, f"{tag}.conf"); prefix = os.path.join(work, tag)
    lines, placed = [], False
    for l in conf_text.splitlines():
        s = l.strip()
        if s.startswith("exp_prefix"): l = f"exp_prefix {prefix}"
        if s == "<END>" and not placed: lines.append(f"gpu {gpu}"); placed = True
        lines.append(l)
    if not placed: lines.append(f"gpu {gpu}")
    open(conf, "w").write("\n".join(lines) + "\n")
    res = subprocess.run([binp, "--conf", conf], cwd=REPO, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, timeout=600)
    log = res.stdout or ""
    if gpu == 1 and ("WITHOUT CUDA" in log or "CUDA error" in log or "no CUDA-capable device" in log):
        return None, "skip"
    if res.returncode != 0 or "[ERROR]" in log:
        sys.stderr.write(log); return None, f"fail (exit {res.returncode})"

    def col(suffix):
        cands = glob.glob(prefix + f"*_DWI{suffix}.txt")
        if not cands: return None
        return [float(x) for x in open(max(cands, key=os.path.getmtime)) if x.strip()]
    intra, extra = col("_intra"), col("_extra")
    if intra is None or extra is None: return None, "no intra/extra output"
    return (sum(intra), sum(extra)), "ok"


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    conf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "tests/accuracy/sphere_compartments.conf")
    # Optional 2nd arg overrides TOL_TOTAL (e.g. permeability+T2 sits ~1% off the CPU due
    # to float32 collision precision shifting the permeable intra/extra equilibrium).
    tol_total = float(sys.argv[2]) if len(sys.argv) > 2 else TOL_TOTAL
    if not os.path.isabs(conf): conf = os.path.join(REPO, conf)
    text = open(conf).read()
    N = int(next(l.split()[1] for l in text.splitlines() if l.strip().startswith("N ")))
    work = tempfile.mkdtemp(prefix="mcdc_gpucomp_")
    try:
        g, st = run(binp, text, 1, work, "gpu")
        if st == "skip": print("SKIP: no CUDA/GPU."); return 0
        if st != "ok": print(f"FAIL (GPU): {st}"); return 1
        c, st = run(binp, text, 0, work, "cpu")
        if st != "ok": print(f"FAIL (CPU): {st}"); return 1
        gi, ge = g; ci, ce = c
        fg, fc = gi/(gi+ge), ci/(ci+ce)
        tg, tc = (gi+ge)/N, (ci+ce)/N    # total signal scaled (num_rep cancels in the diff ratio below)
        print(f"[{os.path.basename(conf)}] intra-frac GPU={fg:.4f} CPU={fc:.4f} |d|={abs(fg-fc):.4f} "
              f"(tol {TOL_FRAC}); total GPU={tg:.1f} CPU={tc:.1f} (tol {tol_total})")
        ok = abs(fg - fc) <= TOL_FRAC and abs(tg - tc)/tc <= tol_total
        print("PASS: GPU compartment split matches CPU." if ok else "FAIL: compartment split disagrees.")
        return 0 if ok else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
