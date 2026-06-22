#!/usr/bin/env python3
"""
GPU volumetric-DWI (voxel subdivisions) validation.

Runs a config with `subdivisions_number`/`subdivisions_file` on GPU (`gpu 1`) and CPU
(`gpu 0`) and checks the per-subdivision signal volume <prefix>_voxels_DWI.txt:
  (a) INVARIANT: for a regular grid the subdivisions partition the voxel, so the GPU's
      per-box signals must sum (over boxes) to the GPU global DWI -- a tight self-check.
  (b) GPU vs CPU: each box's first-acquisition (b~0) signal is ~ its walker population;
      compare the per-box population FRACTION (stable vs independent-RNG placement noise).

Usage: test_gpu_subdivisions.py [config.conf] [tol_frac]
Skips (exit 0) if built without CUDA / no GPU.
"""
import glob, os, shutil, subprocess, sys, tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOL_FRAC = 0.02   # |box-fraction_gpu - box-fraction_cpu|, abs, per box
TOL_SUM  = 1e-4   # relative |sum_boxes - global| on the GPU (partition invariant)


def read_dwi(path):
    return [float(x) for x in open(path) if x.strip()]


def read_voxels(path):
    vals = [float(x) for x in open(path) if x.strip()]
    n_sub, n_rep = int(vals[0]), int(vals[1])
    body = vals[2:]
    return n_sub, n_rep, [body[s*n_rep:(s+1)*n_rep] for s in range(n_sub)]


def sentinel_counts(prefix):
    """CPU sentinel (stuck, illegal). stuck = dropped particles (~0 when healthy);
    illegal = crossing-correction EVENTS (retried, not lost) -- reported, not gated."""
    cands = glob.glob(prefix + "*_simulation_info.txt")
    if not cands: return None
    txt = open(max(cands, key=os.path.getmtime)).read()
    stuck = illegal = 0
    for line in txt.splitlines():
        if "labeled as stuck" in line:            stuck   = int(line.split(":")[-1])
        elif "eliminated due crossings" in line:  illegal = int(line.split(":")[-1])
    return (stuck, illegal)


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
    vox = glob.glob(prefix + "*_voxels_DWI.txt")
    dwi = glob.glob(prefix + "*_DWI.txt")
    dwi = [d for d in dwi if "_voxels" not in d and "_img" not in d]
    if not vox or not dwi: return None, "no voxels/DWI output"
    n_sub, n_rep, boxes = read_voxels(max(vox, key=os.path.getmtime))
    glob_dwi = read_dwi(max(dwi, key=os.path.getmtime))
    return (n_sub, n_rep, boxes, glob_dwi), "ok"


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    conf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "tests/accuracy/subdivisions_free.conf")
    tol_frac = float(sys.argv[2]) if len(sys.argv) > 2 else TOL_FRAC
    if not os.path.isabs(conf): conf = os.path.join(REPO, conf)
    text = open(conf).read()
    N = int(next(l.split()[1] for l in text.splitlines() if l.strip().startswith("N ")))
    work = tempfile.mkdtemp(prefix="mcdc_gpusub_")
    try:
        g, st = run(binp, text, 1, work, "gpu")
        if st == "skip": print("SKIP: no CUDA/GPU."); return 0
        if st != "ok": print(f"FAIL (GPU): {st}"); return 1
        c, st = run(binp, text, 0, work, "cpu")
        if st != "ok": print(f"FAIL (CPU): {st}"); return 1
        (ns_g, nr_g, bx_g, dwi_g) = g
        (ns_c, nr_c, bx_c, dwi_c) = c
        if ns_g != ns_c or nr_g != nr_c:
            print(f"FAIL: shape mismatch GPU({ns_g},{nr_g}) vs CPU({ns_c},{nr_c})"); return 1

        # (a) partition invariant on the GPU: sum_boxes(b0) == global(b0)
        sum_b0 = sum(bx_g[s][0] for s in range(ns_g))
        inv = abs(sum_b0 - dwi_g[0]) / abs(dwi_g[0]) if dwi_g[0] else abs(sum_b0)
        # (b) per-box population fraction GPU vs CPU (use b~0 = first acquisition)
        tot_g = sum(bx_g[s][0] for s in range(ns_g)) or 1.0
        tot_c = sum(bx_c[s][0] for s in range(ns_c)) or 1.0
        worst = 0.0
        for s in range(ns_g):
            fg, fc = bx_g[s][0]/tot_g, bx_c[s][0]/tot_c
            worst = max(worst, abs(fg - fc))
        sc = sentinel_counts(os.path.join(work, "cpu"))
        sc_str = f"stuck={sc[0]} cross-corr-events={sc[1]}" if sc else "n/a"
        print(f"[{os.path.basename(conf)}] n_sub={ns_g} num_rep={nr_g}; "
              f"partition-invariant |sum-global|/global={inv:.2e} (tol {TOL_SUM}); "
              f"worst per-box fraction |d|={worst:.4f} (tol {tol_frac}); CPU {sc_str}")
        ok = inv <= TOL_SUM and worst <= tol_frac
        if sc is not None and sc[0] > 0.01 * N:
            print(f"FAIL: CPU dropped {sc[0]} stuck particles (>1% of N={N}); run is unhealthy.")
            ok = False
        print("PASS: GPU volumetric DWI matches CPU." if ok else "FAIL: volumetric DWI disagrees.")
        return 0 if ok else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
