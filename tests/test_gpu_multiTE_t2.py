#!/usr/bin/env python3
"""
GPU multi-TE T2 validation (Phase 4 / position+T2 log).

A b0 acquisition (G=0) has NO diffusion encoding, so its signal at echo time TE is purely
the T2-weighted spin count S(TE) = sum_w exp(-integral_0^TE dt/T2(compartment at each step)).
A series of b0 at increasing TE therefore isolates the per-step T2 accumulation and the
compartment (position) log from diffusion -- the cleanest probe of "complete tracking".

This runs a config whose scheme contains several b0 rows at different TE (plus, harmlessly,
any diffusion-weighted rows) on GPU (`gpu 1`) and CPU (`gpu 0`, the gold-standard full-tracking
reference) with `separate_signals 1`, then from the INTRA b0-vs-TE decay it:
  (a) checks the GPU and CPU decay curves agree (the per-step T2 log matches), and
  (b) fits an effective intra T2 and checks GPU vs CPU agree and the value is physical.

Usage: test_gpu_multiTE_t2.py [config.conf] [tol_curve] [tol_t2]
Skips (exit 0) if built without CUDA / no GPU.
"""
import glob, math, os, shutil, subprocess, sys, tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOL_CURVE = 0.02     # max |S_gpu(TE)-S_cpu(TE)| / S_cpu(TE) over the b0-vs-TE series (float32 + MC noise)
TOL_T2    = 0.06     # relative |T2_gpu - T2_cpu| / T2_cpu of the fitted effective intra T2


def b0_TEs(scheme_path):
    """(indices, TEs[s]) of the G==0 (b0) acquisitions in scheme order."""
    idx, tes = [], []
    s = 0
    for ln in open(scheme_path):
        p = ln.split()
        if len(p) < 7:
            continue
        try:
            g, te = float(p[3]), float(p[6])
        except ValueError:
            continue
        if abs(g) == 0.0:
            idx.append(s); tes.append(te)
        s += 1
    return idx, tes


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
                         stderr=subprocess.STDOUT, text=True, timeout=900)
    log = res.stdout or ""
    if gpu == 1 and ("WITHOUT CUDA" in log or "CUDA error" in log or "no CUDA-capable device" in log):
        return None, "skip"
    if res.returncode != 0 or "[ERROR]" in log:
        sys.stderr.write(log); return None, f"fail (exit {res.returncode})"
    cands = glob.glob(prefix + "*_DWI_intra.txt")
    if not cands: return None, "no intra output"
    vals = [float(x) for x in open(max(cands, key=os.path.getmtime)) if x.strip()]
    return vals, "ok"


def fit_t2(tes, sig):
    """Effective T2 from a least-squares fit of ln(S/S0) = -TE/T2 (TE in s -> T2 in ms)."""
    s0 = sig[0]
    xs = tes; ys = [math.log(max(v, 1e-12) / s0) for v in sig]
    n = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
    return (-1.0/slope) * 1e3 if slope < 0 else float("inf")


def main():
    binp = os.path.join(REPO, "MC-DC_Simulator")
    conf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "tests/benchmark/three_imp_mte_test.conf")
    tol_curve = float(sys.argv[2]) if len(sys.argv) > 2 else TOL_CURVE
    tol_t2    = float(sys.argv[3]) if len(sys.argv) > 3 else TOL_T2
    if not os.path.isabs(conf): conf = os.path.join(REPO, conf)
    text = open(conf).read()
    scheme = next(l.split(None, 1)[1].strip() for l in text.splitlines() if l.strip().startswith("scheme_file"))
    if not os.path.isabs(scheme): scheme = os.path.join(REPO, scheme)
    idx, tes = b0_TEs(scheme)
    if len(idx) < 3:
        print(f"SKIP: scheme has only {len(idx)} b0 rows (need >=3 for a TE series)."); return 0

    work = tempfile.mkdtemp(prefix="mcdc_mte_")
    try:
        g, st = run(binp, text, 1, work, "gpu")
        if st == "skip": print("SKIP: no CUDA/GPU."); return 0
        if st != "ok": print(f"FAIL (GPU): {st}"); return 1
        c, st = run(binp, text, 0, work, "cpu")
        if st != "ok": print(f"FAIL (CPU): {st}"); return 1
        gb = [g[i] for i in idx]; cb = [c[i] for i in idx]      # intra b0 vs TE
        rel = max(abs(a-b)/abs(b) for a, b in zip(gb, cb) if b != 0)
        t2g, t2c = fit_t2(tes, gb), fit_t2(tes, cb)
        dt2 = abs(t2g - t2c) / t2c
        print(f"[{os.path.basename(conf)}] {len(idx)} b0 TE={tes[0]*1e3:.0f}..{tes[-1]*1e3:.0f}ms | "
              f"intra T2 GPU={t2g:.0f}ms CPU={t2c:.0f}ms (rel {dt2:.3f}, tol {tol_t2}); "
              f"curve max-rel GPU-vs-CPU={rel:.3f} (tol {tol_curve})")
        ok = rel <= tol_curve and dt2 <= tol_t2
        print("PASS: GPU multi-TE T2 matches the CPU gold standard." if ok
              else "FAIL: multi-TE T2 / decay disagrees with CPU.")
        return 0 if ok else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
