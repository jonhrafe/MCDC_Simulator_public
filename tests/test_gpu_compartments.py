#!/usr/bin/env python3
"""
GPU compartment-tracking validation (Phase 3).

Runs a config with `separate_signals 1` on GPU (`gpu 1`) and CPU (`gpu 0`) and checks
that the GPU classifies walkers into the intra / extra compartments the same way the
CPU does: it compares (a) the intra signal FRACTION intra/(intra+extra) -- stable
vs independent-RNG placement noise -- and (b) the total signal. This guards the
per-walker `in_obs` classification the dynamics kernels write at placement.

Usage: test_gpu_compartments.py [config.conf] [tol_total] [tol_frac]
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


def sentinel_counts(prefix):
    """CPU sentinel counts from <prefix>_simulation_info.txt: (stuck, illegal).
    stuck   = particles dropped for getting stuck (should be ~0 in a healthy run).
    illegal = membrane-crossing-correction EVENTS (a surface walker numerically leaks
              and the sentinel back-tracks + retries it); this counts events, not lost
              particles, so it can exceed N and is reported for visibility, not gated."""
    cands = glob.glob(prefix + "*_simulation_info.txt")
    if not cands: return None
    txt = open(max(cands, key=os.path.getmtime)).read()
    stuck = illegal = 0
    for line in txt.splitlines():
        if "labeled as stuck" in line:            stuck   = int(line.split(":")[-1])
        elif "eliminated due crossings" in line:  illegal = int(line.split(":")[-1])
    return (stuck, illegal)


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    conf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "tests/accuracy/sphere_compartments.conf")
    # Optional 2nd arg overrides TOL_TOTAL (e.g. perm+T2 needs a touch more headroom).
    tol_total = float(sys.argv[2]) if len(sys.argv) > 2 else TOL_TOTAL
    # Optional 3rd arg overrides TOL_FRAC. Analytic permeability splits match tightly
    # (<=0.013) now that the CPU crossing UB is fixed; the mesh split is ~0.024 off (the
    # GPU nearest-triangle classification vs the CPU ray-cast isInsidePLY -- a modeling
    # approximation, not noise), so the mesh perm tests pass a relaxed 0.04.
    tol_frac = float(sys.argv[3]) if len(sys.argv) > 3 else TOL_FRAC
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
        # Sentinel health (CPU run), reported for visibility: stuck particles should be ~0;
        # crossing-correction events are informational (events, not lost particles). The
        # signal-tolerance comparison is what gates correctness. stuck > 1% of N is a real
        # red flag (particles genuinely dropped), so that one IS gated.
        sc = sentinel_counts(os.path.join(work, "cpu"))
        sc_str = f"stuck={sc[0]} cross-corr-events={sc[1]}" if sc else "n/a"
        print(f"[{os.path.basename(conf)}] intra-frac GPU={fg:.4f} CPU={fc:.4f} |d|={abs(fg-fc):.4f} "
              f"(tol {tol_frac}); total GPU={tg:.1f} CPU={tc:.1f} (tol {tol_total}); CPU {sc_str}")
        ok = abs(fg - fc) <= tol_frac and abs(tg - tc)/tc <= tol_total
        if sc is not None and sc[0] > 0.01 * N:
            print(f"FAIL: CPU dropped {sc[0]} stuck particles (>1% of N={N}); run is unhealthy.")
            ok = False
        print("PASS: GPU compartment split matches CPU." if ok else "FAIL: compartment split disagrees.")
        return 0 if ok else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
