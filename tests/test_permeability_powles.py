#!/usr/bin/env python3
"""
Powles permeability-model validation.

Checks that the simulator APPLIES the prescribed per-encounter crossing
probabilities correctly: for each permeable obstacle, the empirical
    p_hat = crossings / hits      (per direction)
must match the prescribed `prob_cross_*` within a statistical tolerance.

The experiment (tests/accuracy/permeable_sphere.conf) is a single permeable
sphere with DISTINCT intra/extra diffusivities (Di != De), so the two prescribed
probabilities differ (directional). The test verifies:
  1. p_hat_{i->e} ~= prob_{i->e}  and  p_hat_{e->i} ~= prob_{e->i}  (within K_SIGMA),
  2. the two prescribed probabilities are actually distinct (directional logic on).

This is an IMPLEMENTATION sanity check on easy (smooth-sphere) geometry -- not a
physical-accuracy benchmark of the Powles model on complex meshes.

Counts are read from <prefix>_perm_counters.txt; the run uses num_process 1 so
the (process-shared) counters are exact. Exit 0 = pass, 1 = fail.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

K_SIGMA = 6.0            # tolerance: |p_hat - prob| <= K_SIGMA * sigma
MIN_DIRECTIONAL = 0.2    # require |prob_ie - prob_ei| / max >= this (Di != De)


def read_exp_prefix(text):
    for line in text.splitlines():
        if line.strip().startswith("exp_prefix"):
            return line.split()[1]
    raise ValueError("no exp_prefix line in config")


def write_conf(text, prefix):
    out = []
    for line in text.splitlines():
        out.append(f"exp_prefix {prefix}" if line.strip().startswith("exp_prefix") else line)
    return "\n".join(out) + "\n"


def run_sim(binary, conf_path, workdir):
    res = subprocess.run([binary, "--conf", conf_path], cwd=workdir,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = res.stdout or ""
    if res.returncode != 0 or "[ERROR]" in out or "Assertion" in out:
        sys.stderr.write(out)
        raise RuntimeError(f"simulator failed (exit {res.returncode})")


def parse_counters(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            rows.append(dict(
                type=p[0], id=int(p[1]), kappa=float(p[2]),
                hits_ie=int(p[3]), hits_ei=int(p[4]),
                cross_ie=int(p[5]), cross_ei=int(p[6]),
                prob_ie=float(p[9]), prob_ei=float(p[10])))
    return rows


def check_direction(label, cross, hits, prob):
    if hits == 0:
        return False, f"  FAIL {label}: 0 hits (cannot validate)"
    phat = cross / hits
    sigma = math.sqrt(max(prob * (1.0 - prob), 1e-15) / hits)
    nsig = abs(phat - prob) / sigma if sigma > 0 else float("inf")
    ok = nsig <= K_SIGMA
    return ok, (f"  {'OK  ' if ok else 'FAIL'} {label}: p_hat={phat:.5f} prob={prob:.5f} "
                f"hits={hits} cross={cross} dev={nsig:.2f} sigma")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bin", default=os.path.join(repo_root, "MC-DC_Simulator"))
    p.add_argument("--conf", default=os.path.join(here, "accuracy", "permeable_sphere.conf"))
    p.add_argument("--workdir", default=repo_root)
    args = p.parse_args()

    if not os.path.exists(args.bin):
        print(f"FAIL: simulator binary not found: {args.bin}", file=sys.stderr)
        return 1
    with open(args.conf) as fh:
        conf_text = fh.read()
    base = os.path.basename(read_exp_prefix(conf_text))

    work = tempfile.mkdtemp(prefix="mcdc_perm_")
    try:
        prefix = os.path.join(work, base)
        conf_path = os.path.join(work, "run.conf")
        with open(conf_path, "w") as fh:
            fh.write(write_conf(conf_text, prefix))
        run_sim(args.bin, conf_path, args.workdir)

        counters = prefix + "_perm_counters.txt"
        if not os.path.exists(counters):
            print(f"FAIL: no permeability counters produced ({os.path.basename(counters)}). "
                  f"Is the obstacle permeable and num_process 1?", file=sys.stderr)
            return 1

        rows = parse_counters(counters)
        if not rows:
            print("FAIL: counters file has no obstacle rows", file=sys.stderr)
            return 1

        print(f"Powles permeability check (tolerance {K_SIGMA:g} sigma):")
        ok = True
        for r in rows:
            tag = f"{r['type']}{r['id']}"
            # directional: the two prescribed probabilities must differ (Di != De)
            mx = max(r["prob_ie"], r["prob_ei"], 1e-15)
            rel = abs(r["prob_ie"] - r["prob_ei"]) / mx
            if rel < MIN_DIRECTIONAL:
                ok = False
                print(f"  FAIL {tag}: prob_i_e ({r['prob_ie']:.5f}) and prob_e_i "
                      f"({r['prob_ei']:.5f}) are not distinct (rel {rel:.3f} < {MIN_DIRECTIONAL}); "
                      f"directional logic may be off")
            else:
                print(f"  OK   {tag}: directional prob_i_e={r['prob_ie']:.5f} != "
                      f"prob_e_i={r['prob_ei']:.5f} (rel {rel:.2f})")

            ok_ie, msg_ie = check_direction(f"{tag} i->e", r["cross_ie"], r["hits_ie"], r["prob_ie"])
            ok_ei, msg_ei = check_direction(f"{tag} e->i", r["cross_ei"], r["hits_ei"], r["prob_ei"])
            print(msg_ie); print(msg_ei)
            ok = ok and ok_ie and ok_ei

        if ok:
            print("PASS: empirical crossing probabilities match the prescribed Powles values.")
            return 0
        print("FAIL: empirical crossing probabilities deviate from the prescribed values.")
        return 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
