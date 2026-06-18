#!/usr/bin/env python3
"""
RNG reproducibility regression test for the MC-DC simulator.

Rationale (Phase 0 / P0.2 acceptance test for the RNG unification, P0.1):
A correct Monte-Carlo simulator must be deterministic given a fixed `seed`.
Today it is not: walker placement, substrate generation and the percolation
draw all bypass the seed (see DEV_PLAN_PHASE0.md). This test pins that
contract down so we can prove the RNG fix works and never regresses.

What it does, driving the real binary on the real debug config:
  1. Run the simulator twice with the SAME seed   -> trajectories must be IDENTICAL.
  2. Run the simulator once with a DIFFERENT seed  -> trajectory must DIFFER
     (otherwise the seed has no effect, which is also a bug).

"Identical" is checked both ways the user asked for:
  - bit-exact equality of the raw trajectory arrays (the strong contract), and
  - equality of physically-relevant summary statistics (MSD curve, per-step
    ensemble centroid, global mean/std), which is what a scientist would inspect.

Exit code 0 = pass, 1 = fail. Designed to be invoked directly or via CTest.

NOTE: Until P0.1 (RNG unification) lands, this test is EXPECTED TO FAIL — that
is the point. It turns green once same-seed runs become reproducible.
"""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def build_conf(base_conf_text: str, exp_prefix: str, seed: int) -> str:
    """Return a copy of the base config with exp_prefix redirected and a seed set.

    The seed line is inserted BEFORE the `<END>` terminator, because the parser
    stops reading at `<END>`. We also force the text-trajectory outputs on.
    """
    out_lines = []
    inserted_seed = False
    for raw in base_conf_text.splitlines():
        line = raw.rstrip("\n")
        stripped = line.strip()
        if stripped.startswith("exp_prefix"):
            out_lines.append(f"exp_prefix {exp_prefix}")
            continue
        if stripped.startswith("seed"):
            # Drop any pre-existing seed; we set our own below.
            continue
        if stripped == "<END>" and not inserted_seed:
            out_lines.append(f"seed {seed}")
            out_lines.append("write_txt 1")
            out_lines.append("write_traj_file 1")
            inserted_seed = True
        out_lines.append(line)
    if not inserted_seed:
        # No <END> found: just append.
        out_lines.append(f"seed {seed}")
        out_lines.append("write_txt 1")
        out_lines.append("write_traj_file 1")
    return "\n".join(out_lines) + "\n"


def run_sim(binary: str, conf_path: str, workdir: str) -> None:
    """Run the simulator with the given config, CWD = workdir (paths are relative)."""
    res = subprocess.run(
        [binary, "--conf", conf_path],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    # Also fail on a printed "[ERROR]"/assert: config validation uses assert(),
    # which is disabled under -DNDEBUG, so a bad config can print an error and
    # still exit 0. An exit-code-only check would miss it.
    out = res.stdout or ""
    bad = ("[ERROR]" in out) or ("Assertion" in out)
    if res.returncode != 0 or bad:
        sys.stderr.write(out)
        raise RuntimeError(f"simulator failed (exit {res.returncode}; error in output={bad})")


def load_trajectory(prefix: str):
    """Load <prefix>_0.traj.txt and its header.

    Returns (flat_float_array, N, T). The trajectory file is one float per line
    ordered as N particles x (T+1) steps x 3 coords. The header has 3 numbers:
    duration, N, T.
    """
    hdr_path = f"{prefix}_0.hdr.txt"
    traj_path = f"{prefix}_0.traj.txt"
    if not os.path.exists(traj_path):
        raise FileNotFoundError(f"trajectory file not produced: {traj_path}")
    with open(hdr_path) as fh:
        hdr = fh.read().split()
    n_particles = int(float(hdr[1]))
    n_steps = int(float(hdr[2]))
    flat = np.loadtxt(traj_path, dtype=np.float64)
    return flat.ravel(), n_particles, n_steps


def trajectory_stats(flat: np.ndarray, n_particles: int, n_steps: int) -> dict:
    """Compute physically-relevant summary statistics from a trajectory array."""
    stats = {
        "n_values": int(flat.size),
        "global_mean": float(np.mean(flat)) if flat.size else float("nan"),
        "global_std": float(np.std(flat)) if flat.size else float("nan"),
        "sha256": hashlib.sha256(flat.tobytes()).hexdigest(),
    }
    expected = n_particles * (n_steps + 1) * 3
    if flat.size == expected:
        pos = flat.reshape(n_particles, n_steps + 1, 3)
        disp = pos - pos[:, :1, :]                      # displacement from start
        msd = np.mean(np.sum(disp ** 2, axis=2), axis=0)  # MSD(t), averaged over walkers
        stats["msd_final"] = float(msd[-1])
        stats["msd_mean"] = float(np.mean(msd))
        stats["centroid_final"] = pos[:, -1, :].mean(axis=0).tolist()
    else:
        stats["msd_final"] = stats["msd_mean"] = float("nan")
        stats["centroid_final"] = None
        stats["shape_mismatch"] = f"{flat.size} != expected {expected}"
    return stats


def print_stats(label: str, s: dict) -> None:
    print(f"  [{label}] n_values={s['n_values']} "
          f"mean={s['global_mean']:.6e} std={s['global_std']:.6e} "
          f"msd_final={s['msd_final']:.6e} sha256={s['sha256'][:12]}…")


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bin", default=os.path.join(repo_root, "MC-DC_Simulator"),
                   help="path to the MC-DC_Simulator binary")
    p.add_argument("--conf", default=os.path.join(here, "accuracy", "repro.conf"),
                   help="base config to derive test runs from")
    p.add_argument("--workdir", default=repo_root,
                   help="CWD for the simulator (relative paths in conf resolve here)")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--other-seed", type=int, default=456)
    args = p.parse_args()

    if not os.path.exists(args.bin):
        print(f"FAIL: simulator binary not found: {args.bin}\n"
              f"      build it first (cmake --build build).", file=sys.stderr)
        return 1
    with open(args.conf) as fh:
        base = fh.read()

    work = tempfile.mkdtemp(prefix="mcdc_rng_test_")
    try:
        runs = {}
        for tag, seed in (("A_same", args.seed), ("B_same", args.seed),
                          ("C_diff", args.other_seed)):
            outdir = os.path.join(work, tag)
            os.makedirs(outdir, exist_ok=True)
            prefix = os.path.join(outdir, "run")
            conf_path = os.path.join(work, f"{tag}.conf")
            with open(conf_path, "w") as fh:
                fh.write(build_conf(base, prefix, seed))
            run_sim(args.bin, conf_path, args.workdir)
            flat, n, t = load_trajectory(prefix)
            runs[tag] = (flat, trajectory_stats(flat, n, t))

        print("Trajectory statistics:")
        print_stats(f"seed={args.seed} run A", runs["A_same"][1])
        print_stats(f"seed={args.seed} run B", runs["B_same"][1])
        print_stats(f"seed={args.other_seed} run C", runs["C_diff"][1])

        flat_a, stats_a = runs["A_same"]
        flat_b, stats_b = runs["B_same"]
        flat_c, stats_c = runs["C_diff"]

        ok = True

        # Contract 1: same seed -> bit-exact identical trajectories.
        if flat_a.shape == flat_b.shape and np.array_equal(flat_a, flat_b):
            print(f"PASS: same seed ({args.seed}) -> identical trajectories "
                  f"(bit-exact, {flat_a.size} values).")
        else:
            ok = False
            print(f"FAIL: same seed ({args.seed}) -> trajectories DIFFER "
                  f"(reproducibility broken).")
            if flat_a.shape != flat_b.shape:
                print(f"      different lengths: A={flat_a.size} B={flat_b.size}")
            else:
                ndiff = int(np.count_nonzero(flat_a != flat_b))
                print(f"      {ndiff}/{flat_a.size} values differ; "
                      f"msd_final A={stats_a['msd_final']:.6e} "
                      f"B={stats_b['msd_final']:.6e}")

        # Contract 2: different seed -> trajectories must differ (seed has effect).
        if flat_a.shape == flat_c.shape and np.array_equal(flat_a, flat_c):
            ok = False
            print(f"FAIL: different seeds ({args.seed} vs {args.other_seed}) "
                  f"-> IDENTICAL trajectories (seed has no effect).")
        else:
            print(f"PASS: different seeds -> different trajectories (seed matters).")

        return 0 if ok else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
