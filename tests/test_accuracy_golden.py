#!/usr/bin/env python3
"""
Golden-master accuracy regression test for the MC-DC simulator.

Purpose
-------
Pin down the *core numerical behaviour* of the simulator so that refactors,
optimizations, parallelization changes -- and, later, a GPU re-implementation --
can be checked for "do they still produce the same physics?".

It runs ONE fixed experiment (tests/accuracy/sphere_impermeable.conf):

  * a simple sphere mesh (debug/unitMesh.ply) in a voxel,
  * walkers seeded UNIFORMLY in the voxel -> a natural intra + extra mix,
  * IMPERMEABLE (no permeability) -- the stable invariant that work elsewhere,
    especially the membrane-permeability changes, must never alter,
  * fixed seed and num_process so the run is reproducible.

The diffusion signal it captures -- DWI (real), DWI_intra, DWI_extra -- is the
actual scientific output, so locking it down guards the whole pipeline:
placement, stepping, collision/bouncing, T2, and signal synthesis.

Comparison
----------
Numerical tolerance on each signal:  |run - golden| <= atol + rtol*|golden|.
Defaults are tight (same build, same compiler -> effectively bit-identical).
For a different compiler or a GPU port, loosen with e.g. --rtol 1e-4. The test
reports the worst absolute and relative deviation per signal either way.

Usage
-----
  # check the current build against the committed golden:
  python3 tests/test_accuracy_golden.py
  # (re)generate the golden after an INTENTIONAL, reviewed behaviour change:
  python3 tests/test_accuracy_golden.py --update

Exit code 0 = pass, 1 = fail. Designed to be run directly or via CTest.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

# Signal files the experiment produces (suffixes appended to the exp_prefix).
SIGNALS = ("_DWI.txt", "_DWI_intra.txt", "_DWI_extra.txt")


def read_exp_prefix(conf_text: str) -> str:
    for line in conf_text.splitlines():
        if line.strip().startswith("exp_prefix"):
            return line.split()[1]
    raise ValueError("no exp_prefix line in config")


def write_conf_with_prefix(conf_text: str, new_prefix: str) -> str:
    out = []
    for line in conf_text.splitlines():
        if line.strip().startswith("exp_prefix"):
            out.append(f"exp_prefix {new_prefix}")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def run_sim(binary: str, conf_path: str, workdir: str) -> None:
    res = subprocess.run(
        [binary, "--conf", conf_path],
        cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    # Fail on a non-zero exit OR on any reported error in the output. The latter
    # matters because the simulator's config validation signals failure via
    # assert(), which is a no-op under -DNDEBUG (release builds) and whose return
    # value the caller ignores -- so a broken config can otherwise print "[ERROR]"
    # yet exit 0 and slip past an exit-code-only check.
    out = res.stdout or ""
    bad = ("[ERROR]" in out) or ("Assertion" in out) or ("error:" in out.lower()
                                                          and "no error" not in out.lower())
    if res.returncode != 0 or bad:
        sys.stderr.write(out)
        raise RuntimeError(f"simulator failed (exit {res.returncode}; "
                           f"error in output={bad})")


def load_signal(path: str) -> np.ndarray:
    return np.loadtxt(path, dtype=np.float64).ravel()


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bin", default=os.path.join(repo_root, "MC-DC_Simulator"))
    p.add_argument("--conf", default=os.path.join(here, "accuracy", "sphere_impermeable.conf"))
    p.add_argument("--workdir", default=repo_root,
                   help="CWD for the simulator (relative paths in the conf resolve here)")
    p.add_argument("--golden-dir", default=os.path.join(here, "accuracy", "golden"))
    p.add_argument("--rtol", type=float, default=1e-9,
                   help="relative tolerance (loosen for GPU/other compilers, e.g. 1e-4)")
    p.add_argument("--atol", type=float, default=1e-12, help="absolute tolerance")
    p.add_argument("--update", action="store_true",
                   help="regenerate the golden signals from the current build instead of checking")
    args = p.parse_args()

    if not os.path.exists(args.bin):
        print(f"FAIL: simulator binary not found: {args.bin}\n"
              f"      build it first (cmake --build build).", file=sys.stderr)
        return 1
    with open(args.conf) as fh:
        conf_text = fh.read()
    golden_base = os.path.basename(read_exp_prefix(conf_text))  # e.g. "sphere_impermeable"

    work = tempfile.mkdtemp(prefix="mcdc_accuracy_")
    try:
        prefix = os.path.join(work, golden_base)
        conf_path = os.path.join(work, "run.conf")
        with open(conf_path, "w") as fh:
            fh.write(write_conf_with_prefix(conf_text, prefix))
        run_sim(args.bin, conf_path, args.workdir)

        produced = {sfx: prefix + sfx for sfx in SIGNALS}
        for sfx, path in produced.items():
            if not os.path.exists(path):
                print(f"FAIL: expected output not produced: {os.path.basename(path)}",
                      file=sys.stderr)
                return 1

        if args.update:
            os.makedirs(args.golden_dir, exist_ok=True)
            for sfx, path in produced.items():
                dst = os.path.join(args.golden_dir, golden_base + sfx)
                shutil.copyfile(path, dst)
                print(f"updated golden: {os.path.relpath(dst, repo_root)} "
                      f"({load_signal(path).size} values)")
            print("Golden regenerated. Review the diff before committing.")
            return 0

        ok = True
        print(f"Accuracy check (rtol={args.rtol:g}, atol={args.atol:g}) against "
              f"{os.path.relpath(args.golden_dir, repo_root)}/")
        for sfx, path in produced.items():
            golden_path = os.path.join(args.golden_dir, golden_base + sfx)
            if not os.path.exists(golden_path):
                ok = False
                print(f"  FAIL [{golden_base + sfx}] no golden file; run with --update first.")
                continue
            run = load_signal(path)
            gold = load_signal(golden_path)
            if run.shape != gold.shape:
                ok = False
                print(f"  FAIL [{golden_base + sfx}] shape {run.shape} != golden {gold.shape}")
                continue
            abs_dev = np.abs(run - gold)
            rel_dev = abs_dev / np.maximum(np.abs(gold), np.finfo(float).tiny)
            within = abs_dev <= (args.atol + args.rtol * np.abs(gold))
            n_bad = int(np.count_nonzero(~within))
            tag = "OK  " if n_bad == 0 else "FAIL"
            if n_bad:
                ok = False
            print(f"  {tag} [{golden_base + sfx}] n={run.size} "
                  f"max|abs|={abs_dev.max():.3e} max|rel|={rel_dev.max():.3e} "
                  f"out-of-tol={n_bad}")

        if ok:
            print("PASS: core signal matches golden within tolerance.")
            return 0
        print("FAIL: core signal deviates from golden. If this change is intentional "
              "and reviewed, regenerate with --update.")
        return 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
