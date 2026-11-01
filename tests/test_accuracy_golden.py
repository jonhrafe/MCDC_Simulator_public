#!/usr/bin/env python3
"""
Golden-master accuracy regression test for the MC-DC simulator.

Purpose
-------
Pin down the *core numerical behaviour* of the simulator so that refactors,
optimizations, parallelization changes -- and, later, a GPU re-implementation --
can be checked for "do they still produce the same physics?".

It runs ONE fixed experiment (default: tests/accuracy/sphere_impermeable.conf),
captures the diffusion signal (DWI real / intra / extra -- the actual scientific
output) and compares it to a committed baseline. Locking the signal down guards
the whole pipeline: placement, stepping, collision/bouncing, T2, and synthesis.

Comparison
----------
Numerical tolerance on each signal:  |run - golden| <= atol + rtol*|golden|.
Defaults are tight (same build, same compiler -> effectively bit-identical).
For a different compiler or a GPU port, loosen with e.g. --rtol 1e-4. The test
reports the worst absolute and relative deviation per signal either way.

Usage
-----
  python3 tests/test_accuracy_golden.py            # check against the golden
  python3 tests/test_accuracy_golden.py --update    # regenerate (reviewed change)

Exit code 0 = pass, 1 = fail. Designed to be run directly or via CTest.

Shares its plumbing with the rest of the suite via ``mcdc_testkit`` (config
rewriting, run + error detection, signal loading, tolerance check).
"""

import argparse
import os
import shutil
import sys
import tempfile

# Make ``tests/`` importable when run as a standalone script (CTest passes the
# path explicitly; direct invocation resolves it here).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mcdc_testkit import io, metrics, paths, run  # noqa: E402

# Signal files the experiment produces (suffixes appended to the exp_prefix).
SIGNALS = ("_DWI.txt", "_DWI_intra.txt", "_DWI_extra.txt")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bin", default=None,
                   help="simulator binary (default: $MCDC_BIN or repo-root MC-DC_Simulator)")
    p.add_argument("--conf", default=paths.accuracy("sphere_impermeable.conf"))
    p.add_argument("--workdir", default=paths.REPO_ROOT,
                   help="CWD for the simulator (relative paths in the conf resolve here)")
    p.add_argument("--golden-dir", default=paths.accuracy("golden"))
    p.add_argument("--rtol", type=float, default=1e-9,
                   help="relative tolerance (loosen for GPU/other compilers, e.g. 1e-4)")
    p.add_argument("--atol", type=float, default=1e-12, help="absolute tolerance")
    p.add_argument("--update", action="store_true",
                   help="regenerate the golden signals from the current build instead of checking")
    args = p.parse_args()

    binary = paths.binary_path(args.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}\n"
              f"      build it first (cmake --build build, or ./compile.sh).", file=sys.stderr)
        return 1
    with open(args.conf) as fh:
        conf_text = fh.read()
    golden_base = os.path.basename(io.read_exp_prefix(conf_text))  # e.g. "sphere_impermeable"

    work = tempfile.mkdtemp(prefix="mcdc_accuracy_")
    try:
        prefix = run.run_case(binary, conf_text, work, golden_base, args.workdir)

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
                print(f"updated golden: {os.path.relpath(dst, paths.REPO_ROOT)} "
                      f"({io.load_signal(path).size} values)")
            print("Golden regenerated. Review the diff before committing.")
            return 0

        ok = True
        print(f"Accuracy check (rtol={args.rtol:g}, atol={args.atol:g}) against "
              f"{os.path.relpath(args.golden_dir, paths.REPO_ROOT)}/")
        for sfx, path in produced.items():
            golden_path = os.path.join(args.golden_dir, golden_base + sfx)
            if not os.path.exists(golden_path):
                ok = False
                print(f"  FAIL [{golden_base + sfx}] no golden file; run with --update first.")
                continue
            runv = io.load_signal(path)
            gold = io.load_signal(golden_path)
            if runv.shape != gold.shape:
                ok = False
                print(f"  FAIL [{golden_base + sfx}] shape {runv.shape} != golden {gold.shape}")
                continue
            passed, n_bad, max_abs, max_rel = metrics.within_tol(runv, gold, args.atol, args.rtol)
            if not passed:
                ok = False
            print(f"  {'OK  ' if passed else 'FAIL'} [{golden_base + sfx}] n={runv.size} "
                  f"max|abs|={max_abs:.3e} max|rel|={max_rel:.3e} out-of-tol={n_bad}")

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
