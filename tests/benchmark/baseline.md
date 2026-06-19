# Performance benchmark — baseline (pre-optimization)

A demanding run used to measure the effect of the planned performance work
(Tier 0/1 from the perf review). **Run it with `tests/benchmark/run_bench.sh`**
(from the repo root) before and after each optimization and compare.

## Setup — `porous_media_bench.conf`
- Substrate: `instructions/meshes/PorusMedia.ply` — **39,674 verts / 82,688
  triangles** (`ply_scale 0.001` → a ~13 × 12 × 53 µm porous medium internally).
  Chosen to stress collision detection + the AABB grid.
- Voxel: the mesh bounding box (`file bbox × 1e-6 m`) + ~0.1 µm margin.
- `N = 20000` walkers, `T = 1500` steps, `duration 0.050 s`, `D = 0.6e-9 m²/s`.
- `num_process 10`, `seed 12345`, walkers init **intra** (restricted → heavy
  membrane interaction). Trajectory output OFF; DWI signal ON.

## Baseline numbers — v2.1.0 (commit `ccdb789`), dev host (16 cores, using 10)

| metric | value |
|---|---|
| wall clock | **53.1 s** |
| sim time (avg over 10 workers) | 49 s |
| peak RSS | **122.9 MB** (125,812 kB) |
| DWI `sha256[0:16]` | **`bd4f150c49f4ca98`** |
| DWI b0 (total / intra / extra) | 19829.4 / 19829.4 / 0 |

## After Tier 0 + B4 — v2.1.1.beta (commit `277519a`), same host

| metric | baseline 2.1.0 | + Tier 0 (B1/B2/B3) | + B4 (active-skip) |
|---|---|---|---|
| wall clock | 53.1 s | 49.9 s | **21.2 s** (2.5× vs baseline) |
| peak RSS | 122.9 MB | 132.9 MB | 133.0 MB |
| DWI `sha256[0:16]` | `bd4f150c49f4ca98` | `bd4f150c49f4ca98` | **`bd4f150c49f4ca98`** (bit-exact) |

Notes: B4 (skipping timesteps where the PGSE gradient is off) is the dominant win —
phase/DWI synthesis was ~78% of runtime (measured by A/B: N=12000 30.1 s → 6.7 s
with synthesis disabled). All optimizations are **bit-exact** (identical DWI hash).
RSS is +10 MB from baseline (B1's reused query buffer holds duplicates); the active
list is shared read-only across threads (negligible). Larger peak-RSS reductions
remain Tier-2 work (CSR grid, B5). B2's I/O win is not exercised here (trajectory
output off); it matters when `write_txt`/`write_traj` is on.

## How to interpret before/after
- **Wall time / peak RSS are machine-specific** — a relative reference on the same
  host, not a portable gate.
- **The DWI sha256 is portable and is the correctness reference.** Tier 0 (B1 grid
  scratch buffer, B2 `endl`→`'\n'`, B3 unused VLA) and B6 (per-sim counters) must
  keep it **identical** (`bd4f150c49f4ca98`). B4 (Eigen-batched phase/DWI) may shift
  the last bits via FP reordering — if so, confirm it's noise-level (not a bug),
  then update this hash with a note.
- **Where the wins land:** peak RSS is dominated by the mesh + AABB grid and is
  roughly N-independent, so Tier 0/1 mostly improves **time / allocation churn**;
  peak-RSS reductions come from the deferred CSR-grid / SoA work (B5).
- 10-proc determinism verified (identical DWI hash across repeated runs).
