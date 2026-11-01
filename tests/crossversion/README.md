# Cross-version validation — `tier4_crossversion`

Compare the **current build** against the older **stable** versions (v1.5, v2.0)
and the analytic oracle, per b-shell.

This is the `tier4_crossversion` tier (see [`../README.md`](../README.md#tiers)):
on-demand, mixed hardware, needs the **old binaries**. It answers *"do we still
match the historical stable versions?"*. The **analytic curve is the primary
oracle**; the previous-version signals are a secondary historical overlay.

## Contents

| file | what it is |
|------|-----------|
| `binaries.toml` | registry: `version → {path, argstyle, note, expected_commit}` for every known old build. Parsed with stdlib `tomllib` (python 3.11+). Paths are machine-specific; a missing binary is skipped, never fatal. |
| `../tools/build_old_version.py` | **builds v1.5/v2.0 locally from this repo's own git history** — no external clone. `git worktree add --detach` checks out the release commit (found by bisecting `VERSION_ID` in `src/constants.h`), runs the historical `compile.sh`, caches the binary under `old_binaries/<version>/` (gitignored). `gen_reference.py` calls this automatically the first time a version is needed and just reuses the cache after that (`--force` to rebuild). |
| `../tools/gen_reference.py` | run a chosen build on a per-version config at a fixed seed → freeze `../accuracy/reference_xv/<version>/<case>_DWI.txt` + a `MANIFEST.json` entry. **Refuses** without a fixed positive `--seed` and `--num-process`. |
| `compare_versions.py` | load labelled `_DWI` files + a scheme (+ optional analytic `.npy`) → **markdown report** (per-shell mean S/S0, ADC, NRMSE-vs-reference); **numpy-only, always produced**. PNG plots only if matplotlib imports. |
| `confs/<version>/<case>.conf` | per-version configs (old unit convention) — see [`confs/README.md`](confs/README.md). Only `free_diffusion` is covered today. |

Frozen previous-version signals live under
[`../accuracy/reference_xv/`](../accuracy/reference_xv/) (tiny DWI text files +
`MANIFEST.json`; **never trajectories**).

## Freeze a reference (`gen_reference.py`)

```bash
# Current build on the SI free-diffusion config (the verified mechanism):
python3 tests/tools/gen_reference.py --version current \
    --config tests/accuracy/free_diffusion.conf \
    --seed 12345 --num-process 1 --gpu 0 --case free_diffusion

# v1.5 / v2.0 -- built automatically from this repo's own history on first
# use (tests/tools/build_old_version.py), then cached and just reused:
python3 tests/tools/gen_reference.py --version v1.5 \
    --config tests/crossversion/confs/v1.5/free_diffusion.conf \
    --seed 12345 --num-process 1 --case free_diffusion

python3 tests/tools/gen_reference.py --version v2.0 \
    --config tests/crossversion/confs/v2.0/free_diffusion.conf \
    --seed 12345 --num-process 1 --case free_diffusion
```

- **Determinism guard** — refuses to run unless `--seed` is a fixed **positive**
  integer **and** `--num-process` is given (a MC run is reproducible only for
  both fixed together; default `seed -1` is random).
- **Per-version configs are required for the old builds.** They use different
  `.conf`/unit conventions (`scale_from_stu 1`; voxel bounds always in mm) --
  **do not feed a current SI config to v1.5/v2.0**. Pass the matching config
  under [`confs/<version>/`](confs/README.md); `binaries.toml`'s `argstyle`
  picks the invocation style (both v1.5 and v2.0, built from source here, take
  `--conf`).
- **v1.5/v2.0 binaries build themselves.** No path to hunt down, no external
  clone: the first `gen_reference.py --version v1.5 ...` call builds it from
  this repo's own history and caches it at
  `tests/crossversion/old_binaries/v1.5/MC-DC_Simulator` (gitignored); every
  call after that reuses the cache. Force a rebuild with
  `python3 tests/tools/build_old_version.py --version v1.5 --force`, or build
  ahead of time the same way. There is no v2.2 entry: that line is the current
  code (see [`confs/README.md`](confs/README.md#v22)).

- **Never freezes trajectories** (`write_traj_file` forced to 0; `write_txt`
  forced to 1); only the ~2 KB `_DWI.txt` is copied. If a build wrote only
  `.bfloat`, it is converted to `.txt` on freeze.
- The `MANIFEST.json` entry records: `binary_path`, `git_commit` (via
  `provenance.git_commit` — **null** when the binary is not inside a git repo,
  e.g. v2.0), `version_string` (parsed from `*_simulation_info.txt`; unreliable,
  hence the git commit), `seed`, `num_process`, `N`, `T`, and the
  scheme/conf/signal `sha256`.

## Produce a report (`compare_versions.py`)

```bash
python3 tests/crossversion/compare_versions.py \
    --scheme tests/accuracy/PGSE_sample_scheme.scheme \
    --signal current=tests/accuracy/reference_xv/current/free_diffusion_DWI.txt \
    --signal v1.5=tests/accuracy/reference_xv/v1.5/free_diffusion_DWI.txt \
    --analytic cylinder_r5um_d0p6 \
    --out report.md --plot-dir plots
```

- **Markdown always**, numpy-only: per-shell mean S/S0, per-shell ADC, and
  NRMSE + max|Δ| vs a reference (`--ref`, default `analytic` if present, else the
  first signal). Styled like
  [`../benchmark/three_perm_spheres_report.md`](../benchmark/three_perm_spheres_report.md).
- **N cancels** — S/S0 = signal / mean(signal[b0]), so signals from runs with
  different walker counts are directly comparable (no `--n-walkers`).
- **Plots are optional.** If `matplotlib` is importable (`pip install .[viz]`)
  it writes mean-S/S0-vs-b (linear + log) and ADC-vs-b PNGs; otherwise plots are
  skipped with a note in the report — never a crash.
- **Missing inputs degrade gracefully** — a reference not frozen on this host, a
  wrong-length signal, or an unavailable analytic curve is skipped and noted, not
  fatal.

**How to read the plots**: the analytical reference is
always a solid **black** line; **colour** encodes the simulator version; signals
are never pooled across b-values — one series per shell.

## Determinism caveats

A run is reproducible only for a fixed `seed` **and** fixed `num_process`.
Cross-version and CPU↔GPU results are **never bit-exact** (different RNGs —
`std::mt19937_64` vs cuRAND Philox, and different generators across versions), so
this tier compares **per-shell means within a loose statistical band**, never for
equality.

## CI

- [`../../.github/workflows/crossversion.yml`](../../.github/workflows/crossversion.yml)
  — `workflow_dispatch` only, self-hosted: builds the current binary, freezes it
  with `gen_reference.py`, runs `compare_versions.py`, and uploads the report +
  PNGs as artifacts. v1.5/v2.0 build themselves from this repo's history on
  the runner (needs `git`, `g++`, a working `git worktree`).
- [`../../.github/workflows/nightly.yml`](../../.github/workflows/nightly.yml)
  — the `tier3_analytic` + `gpu` tiers on a self-hosted CUDA runner (not this
  tier, but the sibling on-demand GPU workflow).

## What a real cross-version run needs from you

1. **v1.5 and v2.0 build themselves** — no external binary to source or point
   at. `gen_reference.py --version v1.5|v2.0 ...` builds from this repo's own
   history (`tests/tools/build_old_version.py`) on first use and reuses the
   cached binary after that.
2. A **per-version config** per case (the old unit conventions differ) — see
   [`confs/README.md`](confs/README.md); only `free_diffusion` exists today.
   Freeze each version's `<case>_DWI.txt` with `gen_reference.py`, commit them
   under `../accuracy/reference_xv/<version>/`, then run `compare_versions.py`.
